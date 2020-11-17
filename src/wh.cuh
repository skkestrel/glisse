#pragma once
#include "types.cuh"
#include "data.cuh"
#include "wh.h"

namespace sr
{
namespace wh
{
	using namespace sr::data;

	class WHCudaIntegrator
	{
	public:
		WHIntegrator base;

		Dvf64_3 device_particle_a;
		Dvf64 device_planet_rh;

		uint32_t maxkep;
		float64_t outer_bound;

		Dvf64_3 device_h0_log_0, device_h0_log_1;

		using device_iterator = decltype(device_particle_a.begin());

		void integrate_planets_timeblock(HostPlanetPhaseSpace& pl, float64_t t);
		void integrate_particles_timeblock(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t);
		void gather_particles(const std::vector<size_t>& indices, size_t begin, size_t length);

		void integrate_particles_timeblock_cuda(cudaStream_t stream, size_t planet_data_id, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa);
		void upload_data_cuda(cudaStream_t stream, size_t begin, size_t length);
		void upload_planet_log_cuda(cudaStream_t stream, size_t planet_data_id);
		void recalculate_rh(const HostPlanetPhaseSpace& pl);
		void swap_logs();

		Dvf64_3& device_h0_log(size_t planet_data_id);

		WHCudaIntegrator();
		WHCudaIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config);

		inline device_iterator device_begin()
		{
			return device_particle_a.begin();
		}
	};
}
}
