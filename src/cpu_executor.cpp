#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <thread>

#include "cpu_executor.h"
#include "wh.h"
#include "util.h"
#include "types.h"
#include "convert.h"

namespace sr
{
namespace exec
{
	using namespace sr::data;

	CPUExecutor::CPUExecutor(HostData& _hd, const Configuration& _config, std::ostream& out)
		: hd(_hd), output(out), config(_config) { }

	void CPUExecutor::init()
	{
		sr::convert::to_helio(hd);

		integrator = sr::wh::WHIntegrator(hd.planets, hd.particles, config);
		sr::wh::calculate_planet_metrics(hd.planets, &e_0, nullptr);

		output << std::setprecision(7);
		output << "e_0 (planets) = " << e_0 << std::endl;
		output << "n_particle = " << hd.particles.n() << std::endl;
		output << "n_particle_alive = " << hd.particles.n_alive() << std::endl;
		output << "==================================" << std::endl;
		output << "Sending initial conditions to GPU." << std::endl;

		starttime = std::chrono::high_resolution_clock::now();
		output << "       Starting simulation.       " << std::endl << std::endl;

		step_planets();
		swap_logs();
	}

	void CPUExecutor::swap_logs()
	{
		hd.planets.swap_logs();
		integrator.swap_logs();
	}

	void CPUExecutor::step_planets()
	{
		integrator.integrate_planets_timeblock(hd.planets, t);
	}

	void CPUExecutor::add_job(const std::function<void()>& job)
	{
		work.push_back(std::move(job));
	}

	double CPUExecutor::time() const
	{
		auto now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> millis = now - starttime;
		return millis.count() / 60000;
	}

	void CPUExecutor::loop(double* cputimeout)
	{
		for (auto& i : work) i();
		work.clear();

		std::vector<std::thread> threads;
		if (config.num_thread > 1)
		{
			for (size_t i = 0; i < config.num_thread; i++)
			{
				threads.push_back(std::thread([i, this]()
					{
						size_t total = hd.particles.n_alive();
						integrator.integrate_particles_timeblock(
								hd.planets,
								hd.particles,
								total * i / config.num_thread,
								// The length is not exactly total / config.num_thread due to roundoff
								total * (i + 1) / config.num_thread - total * i / config.num_thread,
								t);
					}));
			}
		}
		else
		{
			size_t total = hd.particles.n_alive();
			integrator.integrate_particles_timeblock(
					hd.planets,
					hd.particles,
					0,
					total,
					t);
		}

		// Copy ctor
		hd.planets_snapshot = hd.planets.base;

		t += config.dt * static_cast<double>(config.tbsize);

		step_planets();
		for (auto& i : threads)
		{
			i.join();
		}
		swap_logs();

		resync();

#pragma GCC warning "TODO"
		*cputimeout = 0;
	}

	void CPUExecutor::resync()
	{
		if (hd.particles.n_alive() == 0)
		{
			return;
		}

		size_t prev_alive = hd.particles.n_alive();

		auto gather_indices = hd.particles.stable_partition_unflagged(0, prev_alive);
		integrator.gather_particles(*gather_indices, 0, prev_alive);

		size_t diff = prev_alive - hd.particles.n_alive();

		for (size_t i = hd.particles.n_alive(); i < hd.particles.n_alive() + diff; i++)
		{
			if ((hd.particles.deathflags()[i] & 0x00FF) == 0x0001)
			{
				hd.particles.deathflags()[i] |= 0x0080;
			}
		}

		gather_indices = hd.particles.stable_partition_alive(prev_alive - diff, diff);
		integrator.gather_particles(*gather_indices, prev_alive - diff, diff);
	}

	void CPUExecutor::finish()
	{
		for (auto& i : work) i();
		work.clear();

		resync();

		for (auto& i : work) i();
		work.clear();

		output << "Simulation finished. t = " << t << ". n_particle = " << hd.particles.n_alive() << std::endl;
	}
}
}
