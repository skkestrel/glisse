#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <thread>

#include "util.cuh"
#include "util.h"
#include "types.h"
#include "executor.cuh"
#include "wh.cuh"
#include "convert.h"

namespace sr
{
namespace exec
{
	using namespace sr::wh;
	using namespace sr::util;
	using namespace sr::convert;
	using namespace sr::data;

	ExecutorData::ExecutorData() { }
	ExecutorData::ExecutorData(size_t n)
	{
		r = v = std::vector<f64_3>(n);
		deathflags = std::vector<uint16_t>(n);
		id = std::vector<uint32_t>(n);
		deathtime_index = std::vector<uint32_t>(n);
	}

	struct DeviceParticleUnflaggedPredicate
	{
		template<typename Tuple>
		__host__ __device__
		bool operator()(const Tuple& args)
		{
			uint16_t flag = thrust::get<2>(thrust::get<0>(args));
			return flag == 0;
		}
	};

	Executor::Executor(HostData& _hd, DeviceData& _dd, const Configuration& _config, std::ostream& out)
		: hd(_hd), dd(_dd), output(out), config(_config) { }

	void Executor::init()
	{
		if (!config.use_gpu)
		{
			output << "Executable was compiled with CUDA but USE_GPU was disabled!" << std::endl;
			throw std::exception();
		}

		to_helio(hd);

		integrator = sr::wh::WHCudaIntegrator(hd.planets, hd.particles, config);
		calculate_planet_metrics(hd.planets, &e_0, nullptr);

		output << std::setprecision(7);
		output << "e_0 (planets) = " << e_0 << std::endl;
		output << "n_particle = " << hd.particles.n() << std::endl;
		output << "n_particle_alive = " << hd.particles.n_alive() << std::endl;
		output << "==================================" << std::endl;
		output << "Sending initial conditions to GPU." << std::endl;

		cudaStreamCreate(&main_stream);
		cudaStreamCreate(&htd_stream);
		cudaStreamCreate(&par_stream);
		cudaStreamCreate(&dth_stream);

		cudaEventCreate(&start_event);
		cudaEventCreate(&cpu_finish_event);
		cudaEventCreate(&gpu_finish_event);

		dd.particles = DeviceParticlePhaseSpace(hd.particles.n());

		dd.planets0 = DevicePlanetPhaseSpace(hd.planets.n(), config.tbsize);
		dd.planets1 = DevicePlanetPhaseSpace(hd.planets.n(), config.tbsize);
		dd.planet_data_id = 0;

		memcpy_htd(dd.planet_phase_space().m, hd.planets.m(), htd_stream);
		cudaStreamSynchronize(htd_stream);
		dd.planet_data_id++;
		memcpy_htd(dd.planet_phase_space().m, hd.planets.m(), htd_stream);
		cudaStreamSynchronize(htd_stream);

		if (hd.particles.n() > 0)
		{
			upload_data(0, hd.particles.n());
		}

		download_data();

		starttime = std::chrono::high_resolution_clock::now();

		output << "       Starting simulation.       " << std::endl << std::endl;

		if (encounter_output)
		{
			*encounter_output << std::setprecision(17);
		}

		step_and_upload_planets();
	}

	void Executor::swap_logs()
	{
		hd.planets.swap_logs();
		integrator.swap_logs();
	}

	void Executor::step_and_upload_planets()
	{
		integrator.integrate_planets_timeblock(hd.planets, t);

		swap_logs();

		// We only upload the planet log if any particles are going to use the planet log on the GPU
		// Cases where the planet log is not used by the particles:
		// - There are no particles alive on the GPU, AND there are no particles in close encounters on the CPU
		// since the particles that survive close encounters can make it to the GPU at the end of this timestep
		// and thus the next planet chunk will be required

		if (dd.particle_phase_space().n_alive > 0 || hd.particles.n_encounter() > 0)
		{
			upload_planet_log();
		}
	}

	void Executor::upload_data(size_t begin, size_t length)
	{
		auto& particles = dd.particle_phase_space();
		particles.n_alive = hd.particles.n_alive();
		integrator.upload_data_cuda(htd_stream, begin, length);

		memcpy_htd(particles.r, hd.particles.r(), htd_stream, begin, begin, length);
		cudaStreamSynchronize(htd_stream);
		memcpy_htd(particles.v, hd.particles.v(), htd_stream, begin, begin, length);
		cudaStreamSynchronize(htd_stream);
		memcpy_htd(particles.deathflags, hd.particles.deathflags(), htd_stream, begin, begin, length);
		cudaStreamSynchronize(htd_stream);
		memcpy_htd(particles.id, hd.particles.id(), htd_stream, begin, begin, length);
		cudaStreamSynchronize(htd_stream);
	}

	void Executor::add_job(const std::function<void()>& job)
	{
		work.push_back(std::move(job));
	}

	void Executor::download_data(bool ignore_errors)
	{
		auto& particles = dd.particle_phase_space();

		Vu32 prev_ids(hd.particles.id().begin(), hd.particles.id().end());

		memcpy_dth(hd.particles.r(), particles.r, dth_stream, 0, 0, particles.n_alive);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(hd.particles.v(), particles.v, dth_stream, 0, 0, particles.n_alive);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(hd.particles.id(), particles.id, dth_stream, 0, 0, particles.n_alive);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(hd.particles.deathflags(), particles.deathflags, dth_stream, 0, 0, particles.n_alive);
		cudaStreamSynchronize(dth_stream);

		// This should NEVER happen. I think this is a recoverable 
		// error, by swapping particle indices on the host, but that sounds annoying...
		if (prev_ids != hd.particles.id())
		{
			output << "WARNING! ID MISMATCH! WARNING!" << std::endl;

			if (!ignore_errors)
			{
				throw std::exception();
			}
		}

		hd.particles.n_alive() = dd.particle_phase_space().n_alive;
	}

	void Executor::upload_planet_log()
	{
		dd.planet_data_id++;
		auto& planets = dd.planet_phase_space();

		memcpy_htd(planets.r_log, hd.planets.r_log().slow, htd_stream);
		cudaStreamSynchronize(htd_stream);

		integrator.upload_planet_log_cuda(htd_stream, dd.planet_data_id);
	}


	double Executor::time() const
	{
		auto now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> millis = now - starttime;
		return millis.count() / 60000;
	}

	void Executor::loop(double* cputimeout, double* gputimeout)
	{
		std::thread cpu_thread;
		
		if (dd.particle_phase_space().n_alive > 0)
		{
			cudaEventRecord(start_event, main_stream);
			integrator.integrate_particles_timeblock_cuda(main_stream, dd.planet_data_id, dd.planet_phase_space(), dd.particle_phase_space());
			cudaEventRecord(gpu_finish_event, main_stream);
		}

		// The queued work should begin RIGHT after the CUDA call
		for (auto& i : work) i();
		work.clear();


		size_t encounter_start = hd.particles.n_alive() - hd.particles.n_encounter();
		for (size_t i = encounter_start; i < hd.particles.n_alive(); i++)
		{
			integrator.integrate_encounter_particle_catchup(hd.planets, hd.particles, i,
					ed.deathtime_index[i - encounter_start],
					t - config.dt * static_cast<double>(config.tbsize - ed.deathtime_index[i - encounter_start])
				);
		}

		auto gather_indices = hd.particles.stable_partition_alive(encounter_start, hd.particles.n_encounter());
		integrator.gather_particles(*gather_indices, encounter_start, hd.particles.n_encounter());
		upload_data(encounter_start, hd.particles.n_encounter());

		// Fill deathtime index with 0 so that the continuation will work
		thrust::fill(thrust::cuda::par.on(htd_stream), dd.particles.deathtime_index.begin() + encounter_start,
				dd.particles.deathtime_index.begin() + encounter_start + hd.particles.n_encounter(), 0);

		// The snapshot contains the planet states at the end of the previous timestep - 
		// consider removing this? We can use hd.planets.*_log_old()[-1] to replicate this functionality

		// Copy assignment ctor
		hd.planets_snapshot = hd.planets.base;

		t += config.dt * static_cast<double>(config.tbsize);
		step_and_upload_planets();

		if (dd.particle_phase_space().n_alive > 0)
		{
			cudaStreamSynchronize(htd_stream);
			cudaEventRecord(cpu_finish_event, par_stream);
			cudaEventSynchronize(gpu_finish_event);

			float cputime, gputime;
			cudaEventElapsedTime(&cputime, start_event, cpu_finish_event);
			cudaEventElapsedTime(&gputime, start_event, gpu_finish_event);
			if (cputimeout) *cputimeout = cputime;
			if (gputimeout) *gputimeout = gputime;

			// There's nothing to resync if all the particles on the device are dead!
			// Although dd.particles.n_alive can be out of sync with dd.particles.deathflags before
			// resync() is called, this is safe:
			// - The MVS kernel does not revive particles, so resync() will never INCREASE n_alive
			// - dd.particles.n_alive is adjusted by the close encounter handler just BEFORE this call
			resync();
		}
	}

	void Executor::resync()
	{
		auto& particles = dd.particle_phase_space();
		size_t prev_alive = particles.n_alive;

		auto partition_it = thrust::make_zip_iterator(thrust::make_tuple(particles.begin(), integrator.device_begin()));
		particles.n_alive = thrust::stable_partition(thrust::cuda::par.on(main_stream),
				partition_it, partition_it + particles.n_alive, DeviceParticleUnflaggedPredicate()) - partition_it;
		cudaStreamSynchronize(main_stream);

		size_t diff = prev_alive - particles.n_alive;

		ed = ExecutorData(diff);

		memcpy_dth(ed.r, particles.r, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(ed.v, particles.v, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(ed.id, particles.id, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(ed.deathtime_index, particles.deathtime_index, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(ed.deathflags, particles.deathflags, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);

		for (size_t i = 0; i < diff; i++)
		{
			if ((ed.deathflags[i] & 0x00FF) == 0x0001)
			{
				if (config.resolve_encounters)
				{
				}
				else
				{
					// If encounters are not being dealt with, kill the particle!
					ed.deathflags[i] |= 0x0080;
				}
			}
		}

		std::unique_ptr<std::vector<size_t>> ed_indices;
		stable_partition_alive_indices(ed.deathflags, 0, diff, &ed_indices);
		gather(ed.r, *ed_indices, 0, diff);
		gather(ed.v, *ed_indices, 0, diff);
		gather(ed.id, *ed_indices, 0, diff);
		gather(ed.deathflags, *ed_indices, 0, diff);
		gather(ed.deathtime_index, *ed_indices, 0, diff);

		std::unordered_map<size_t, size_t> indices;
		for (size_t i = 0; i < prev_alive; i++)
		{
			indices[hd.particles.id()[i]] = i;
		}

		for (size_t i = 0; i < diff; i++)
		{
			size_t index = indices[ed.id[i]];
			hd.particles.r()[index] = ed.r[i];
			hd.particles.v()[index] = ed.v[i];
			hd.particles.deathflags()[index] = ed.deathflags[i];

			if (ed.deathflags[i])
			{
				hd.particles.deathtime()[index] = static_cast<float>(t - config.dt * static_cast<double>(config.tbsize - ed.deathtime_index[i]));
			}
		}

		auto gather_indices = hd.particles.stable_partition_alive(0, prev_alive);
		integrator.gather_particles(*gather_indices, 0, prev_alive);

		hd.particles.n_encounter() = hd.particles.n_alive() - particles.n_alive;

		size_t encounter_start = particles.n_alive;

		add_job([encounter_start, diff, this]()
			{
				if (encounter_output)
				{
					for (size_t i = hd.particles.n_encounter(); i < diff; i++)
					{
						*encounter_output << hd.particles.r()[encounter_start + i] << std::endl;
						*encounter_output << hd.particles.v()[encounter_start + i] << std::endl;
						*encounter_output << hd.particles.id()[encounter_start + i] << " "
							<< hd.particles.deathflags()[encounter_start + i] << " "
							<< t - config.dt * static_cast<double>(config.tbsize - ed.deathtime_index[i]) << " death"
							<< std::endl;
						*encounter_output << hd.planets.n_alive() << std::endl;

						*encounter_output << hd.planets.m()[0] << std::endl;
						*encounter_output << f64_3(0) << std::endl;
						*encounter_output << f64_3(0) << std::endl;
						*encounter_output << hd.planets.id()[0] << std::endl;
						for (size_t j = 1; j < hd.planets.n_alive(); j++)
						{
							*encounter_output << hd.planets.m()[j] << std::endl;
							*encounter_output << hd.planets.r_log().slow[ed.deathtime_index[i] * (hd.planets.n() - 1) + j - 1] << std::endl;
							*encounter_output << hd.planets.v_log().slow[ed.deathtime_index[i] * (hd.planets.n() - 1) + j - 1] << std::endl;
							*encounter_output << hd.planets.id()[j] << std::endl;
						}
					}

					*encounter_output << std::flush;
				}
			});
	}


	void Executor::finish()
	{
		cudaStreamSynchronize(main_stream);

		for (auto& i : work) i();
		work.clear();

		resync();

		for (auto& i : work) i();
		work.clear();

		output << "Simulation finished. t = " << t << ". n_particle = " << hd.particles.n_alive() << std::endl;
	}
}
}
