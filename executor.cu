#include <iomanip>
#include <fstream>
#include <algorithm>

#include "executor.h"
#include "convert.h"
#include "wh.h"

struct DeviceParticleUnflaggedPredicate
{
	template<typename Tuple>
	__host__ __device__
	bool operator()(const Tuple& args)
	{
		uint8_t flag = thrust::get<3>(args);
		return flag == 0;
	}
};

Executor::Executor(HostData& hd, DeviceData& dd, std::ostream& out) : hd(hd), dd(dd), print_every(10), print_counter(0), tbsize(128), output(out), timing_output(nullptr) { }

void Executor::init()
{
	to_helio(hd);

	calculate_planet_metrics(hd.planets, &e_0, nullptr);

	output << std::setprecision(17);
	output << "e_0 (planets) = " << e_0 << std::endl;
	output << "t_0 = " << t << std::endl;
	output << "dt = " << dt << std::endl;
	output << "t_f = " << t_f << std::endl;
	output << "n_particle = " << hd.particles.n << std::endl;
	output << "==================================" << std::endl;
	output << "Running for half a time step.     " << std::endl;

	initialize(hd.planets, hd.particles);

	output << "==================================" << std::endl;
	output << "Sending initial conditions to GPU." << std::endl;

	cudaStreamCreate(&main_stream);
	cudaStreamCreate(&htd_stream);
	cudaStreamCreate(&par_stream);
	cudaStreamCreate(&dth_stream);

	upload_data();
	cudaStreamSynchronize(htd_stream);

	download_data();
	cudaStreamSynchronize(dth_stream);

	starttime = std::chrono::high_resolution_clock::now();

	output << "       Starting simulation.       " << std::endl << std::endl;

	if (timing_output)
	{
		*timing_output << std::setprecision(17);
	}
}

void Executor::upload_data()
{
	dd.particles0 = DeviceParticlePhaseSpace(hd.particles.n);
	dd.particles1 = DeviceParticlePhaseSpace(hd.particles.n);

	dd.planets0 = DevicePlanetPhaseSpace(hd.planets.n, hd.tbsize);
	dd.planets1 = DevicePlanetPhaseSpace(hd.planets.n, hd.tbsize);

	dd.planet_data_id = 0;
	dd.particle_data_id = 0;

	auto& particles = dd.particle_phase_space();
	thrust::copy(thrust::cuda::par.on(htd_stream), hd.particles.r.begin(), hd.particles.r.end(), particles.r.begin());
	thrust::copy(thrust::cuda::par.on(htd_stream), hd.particles.v.begin(), hd.particles.v.end(), particles.v.begin());
	thrust::copy(thrust::cuda::par.on(htd_stream), hd.particles.a.begin(), hd.particles.a.end(), particles.a.begin());
	thrust::copy(thrust::cuda::par.on(htd_stream), hd.particles.id.begin(), hd.particles.id.end(), particles.id.begin());

	thrust::copy(thrust::cuda::par.on(htd_stream), hd.planets.m.begin(), hd.planets.m.end(), dd.planet_phase_space().m.begin());
	dd.planet_data_id++;
	thrust::copy(thrust::cuda::par.on(htd_stream), hd.planets.m.begin(), hd.planets.m.end(), dd.planet_phase_space().m.begin());
	
	cudaStreamSynchronize(htd_stream);
}

void Executor::download_data()
{
	auto& particles = dd.particle_phase_space();
	size_t n = hd.particles.n;

	Hvu8 flags(n);

	thrust::copy(thrust::cuda::par.on(dth_stream), particles.r.begin(), particles.r.begin() + n, hd.particles.r.begin());
	thrust::copy(thrust::cuda::par.on(dth_stream), particles.v.begin(), particles.v.begin() + n, hd.particles.v.begin());
	thrust::copy(thrust::cuda::par.on(dth_stream), particles.a.begin(), particles.a.begin() + n, hd.particles.a.begin());
	thrust::copy(thrust::cuda::par.on(dth_stream), particles.flags.begin(), particles.flags.begin() + n, flags.begin());
	thrust::copy(thrust::cuda::par.on(dth_stream), particles.id.begin(), particles.id.begin() + n, hd.particles.id.begin());
	cudaStreamSynchronize(dth_stream);


	for (size_t i = 0; i < n; i++)
	{
		hd.particles.flags[i] = flags[i];
	}
	hd.particles.n_alive = dd.particle_phase_space().n_alive;

	/*
	// zip will crash the program

	auto iterator = thrust::make_zip_iterator(thrust::make_tuple(
				ps.r.begin(),
				ps.v.begin(),
				ps.flags.begin(), ps.deathtime.begin(), ps.id.begin()));
	thrust::copy(thrust::cuda::par.on(stream),
			iterator,
			iterator + n,
			thrust::make_zip_iterator(thrust::make_tuple(
					hd.particles.r.begin(), hd.particles.v.begin(),
					hd.particles.flags.begin(), hd.deathtime.begin(), hd.id.begin())));
	 */
}

void Executor::upload_planet_log()
{
	dd.planet_data_id++;
	auto& planets = dd.planet_phase_space();
	thrust::copy(thrust::cuda::par.on(htd_stream), hd.planets.h0_log.begin(), hd.planets.h0_log.end(), planets.h0_log.begin());
	thrust::copy(thrust::cuda::par.on(htd_stream), hd.planets.r_log.begin(), hd.planets.r_log.end(), planets.r_log.begin());
}

double Executor::time() const
{
	auto now = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> millis = now - starttime;
	return millis.count() / 60000;
}

void Executor::loop()
{
	if (print_counter % print_every == 0)
	{
		double e;
		calculate_planet_metrics(hd.planets, &e, nullptr);

		double elapsed = time();
		double total = elapsed * (t_f - t_0) / (t - t_0);

		output << "t=" << t << " (" << elapsed / total * 100 << "% " << elapsed << "m elapsed, " << total << "m total " << total - elapsed << "m remain)" << std::endl;
		output << "Error = " << (e - e_0) / e_0 * 100 << ", " << dd.particle_phase_space().n_alive << " particles remaining" << std::endl;
	}
	print_counter++;

	// advance planets
	for (size_t i = 0; i < tbsize; i++)
	{
		step_planets(hd.planets, t, i, dt);

		if (timing_output)
		{
			double e_;
			f64_3 l_;
			calculate_planet_metrics(hd.planets, &e_, &l_);
		
			*timing_output << "ep " << e_ << std::endl;
			*timing_output << "lp " << l_.x << " " << l_.y << " " << l_.z << std::endl;
		}
	}

	upload_planet_log();
	cudaStreamSynchronize(htd_stream);

	download_data();
	cudaStreamSynchronize(dth_stream);

	cudaStreamSynchronize(main_stream);
	resync();

	step_particles_cuda(main_stream, dd.planet_phase_space(), dd.particle_phase_space(), tbsize, dt);

	t += dt * tbsize;
}

std::ofstream anilog("animation.out");
void Executor::animate()
{
	Vf64_3 r(hd.planets.n);
	Vf64_3 v(hd.planets.n);
	Vf64_3 rp(hd.particles.n);
	Vf64_3 vp(hd.particles.n);

	f64_3 bary_r, bary_v;

	find_barycenter(hd.planets.r, hd.planets.v, hd.planets.m, hd.planets.n, bary_r, bary_v);

	for (size_t i = 0; i < hd.planets.n; i++)
	{
		r[i] = hd.planets.r[i] - bary_r;
		v[i] = hd.planets.v[i] - bary_v;
	}
	for (size_t i = 0; i < hd.particles.n; i++)
	{
		rp[i] = hd.particles.r[i] - bary_r;
		vp[i] = hd.particles.v[i] - bary_v;
	}

	for (size_t j = 0; j < hd.planets.n; j++)
	{
		anilog << std::setprecision(17) << r[j].x << " " << r[j].y << " " << r[j].z << std::endl << std::flush;
		anilog << std::setprecision(17) << v[j].x << " " << v[j].y << " " << v[j].z << std::endl << std::flush;
	}
/*
	for (size_t j = 0; j < hd.particles.n; j++)
	{
		anilog << rp[j].x << " " << rp[j].y << " " << rp[j].z << std::endl << std::flush;
	}
*/
}

void Executor::resync()
{
	auto& particles = dd.particle_phase_space();
	particles.n_alive = thrust::partition(thrust::cuda::par.on(par_stream), particles.begin(), particles.begin() + particles.n_alive, DeviceParticleUnflaggedPredicate())
		- particles.begin();

	// partition data on GPU
	// start GPU kernel
	// pull data to CPU
	// run CEs
	// reorder: alive particles (just finished CEs) first, dead particles last
	// push data to last bit in GPU
	// edit gpu n_alive
	// repeat
	
	/*
	DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;
	DevicePhaseSpace& other_ps = dd.phase_space_id % 2 ? dd.ps1 : dd.ps0;

	thrust::copy(thrust::cuda::par.on(main_stream), ps.flags.begin(), ps.flags.begin() + hd.particles.n_alive, hd.part_flags.begin());
	Hvu32 indices = Hvu32(hd.n_part);

	// start placing alive particles at the front and dead particles at the back
	uint32_t front_counter = 0;
	uint32_t back_counter = hd.particles.n_alive - 1;

	for (size_t i = 0; i < hd.particles.n_alive; i++)
	{
		if (hd.part_flags[i] & 0x0001) // dead particle
		{
			indices[back_counter--] = i;
		}
		else
		{
			indices[front_counter++] = i;
		}
	}
	for (size_t i = hd.particles.n_alive; i < hd.n_part; i++)
	{
		indices[i] = i;
	}

	cudaStreamSynchronize(main_stream);
	thrust::copy(thrust::cuda::par.on(main_stream), indices.begin(), indices.end(), dd.gather_indices.begin());

	size_t resyncd = hd.particles.n_alive - front_counter;
	hd.particles.n_alive = front_counter;

	cudaStreamSynchronize(main_stream);
	cudaStreamSynchronize(main_stream);

	thrust::gather(thrust::cuda::par.on(main_stream),
			dd.gather_indices.begin(),
			dd.gather_indices.end(),
			thrust::make_zip_iterator(thrust::make_tuple(
					ps.r.begin(), ps.v.begin(),
					ps.flags.begin(), ps.deathtime.begin(), ps.id.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(
					other_ps.r.begin(), other_ps.v.begin(),
					other_ps.flags.begin(), other_ps.deathtime.begin(), other_ps.id.begin())));

	dd.phase_space_id++;
*/
}


void Executor::finish()
{
	cudaStreamSynchronize(main_stream);
	resync();

	cudaStreamSynchronize(par_stream);

	output << "Simulation finished. t = " << t << ". n_particle = " << hd.particles.n_alive << std::endl;

	cudaStreamSynchronize(dth_stream);
	download_data();
}
