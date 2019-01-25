#pragma once
#include "types.cuh"

namespace sr
{
namespace data
{
	struct DeviceParticlePhaseSpace
	{
		Dvf64_3 r, v;

		Dvu16 deathflags;
		Dvu32 deathtime_index;
		Dvu32 id;

		Dvu32 gather_indices;

		size_t n_total;
		size_t n_alive;

		DeviceParticlePhaseSpace(const DeviceParticlePhaseSpace&) = delete;

		inline DeviceParticlePhaseSpace() { }
		inline DeviceParticlePhaseSpace(size_t n) : r(n), v(n), deathflags(n), deathtime_index(n), id(n), n_total(n), n_alive(n), gather_indices(n) { }

		using iterator_tuple = decltype(thrust::make_tuple(r.begin(), v.begin(), deathflags.begin(), deathtime_index.begin(), id.begin()));
		using iterator = thrust::zip_iterator<iterator_tuple>;

		inline iterator begin()
		{
			return thrust::make_zip_iterator(thrust::make_tuple(
				r.begin(), v.begin(), deathflags.begin(), deathtime_index.begin(), id.begin()));
		}
	};

	struct DevicePlanetPhaseSpace
	{
		Dvf64 m;
		Dvf64_3 r_log;

		size_t n_total;
		size_t n_alive;

		DevicePlanetPhaseSpace(const DevicePlanetPhaseSpace&) = delete;

		inline DevicePlanetPhaseSpace() { }
		inline DevicePlanetPhaseSpace(size_t n, size_t tbsize) : m(n), r_log(n * tbsize), n_total(n), n_alive(n) { }
	};

	struct DeviceData
	{
		// double buffer for HtD transfer of planet locations
		DeviceParticlePhaseSpace particles;
		DevicePlanetPhaseSpace planets0, planets1;

		size_t planet_data_id;

		DeviceData(const DeviceData&) = delete;
		inline DeviceData() { }
		inline DeviceParticlePhaseSpace& particle_phase_space() { return particles; }
		inline DevicePlanetPhaseSpace& planet_phase_space() { return planet_data_id % 2 ? planets1 : planets0; }
	};
}
}
