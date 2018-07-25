#pragma once
#include "integrator.h"

namespace sr
{
namespace integrator
{
	using namespace sr::data;

	class CudaIntegrator : public Integrator
	{
	public:
		virtual void integrate_particles_timeblock_cuda(cudaStream_t stream, const HostPlanetPhaseSpace& pl_h, size_t planet_data_id, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa) = 0;
		virtual void upload_planet_log_cuda(cudaStream_t stream, size_t planet_data_id) = 0;
		virtual void upload_data_cuda(cudaStream_t stream, size_t begin, size_t length) = 0;
	};
}
}
