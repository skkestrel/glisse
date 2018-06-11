#pragma once
#include "types.h"
#include <string>

#include <cstdlib>
using size_t = std::size_t;

struct DeviceParticlePhaseSpace
{
	Dvf64_3 r, v;
	Dvu8 flags;

	size_t n_total;
	size_t n_alive;
	size_t n_encounter;

	inline DeviceParticlePhaseSpace() { }
	inline DeviceParticlePhaseSpace(size_t n) : r(n), v(n), flags(n), n_total(n), n_alive(n) { }
};

struct DevicePlanetPhaseSpace
{
	Dvf64 m;
	Dvf64_3 r_log;
	Dvf64_3 h0_log;

	size_t n_total;
	size_t n_alive;

	inline DevicePlanetPhaseSpace() { }
	inline DevicePlanetPhaseSpace(size_t n, size_t tbsize) : m(n), r_log(n * tbsize), h0_log(tbsize),
		n_total(n), n_alive(n) { }
};

struct DeviceData
{
	// double buffer for HtD transfer of planet locations
	DeviceParticlePhaseSpace particles0, particles1;
	DevicePlanetPhaseSpace planets0, planets1;

	size_t particle_data_id, planet_data_id;

	inline DeviceData() { }
	inline DeviceParticlePhaseSpace& particle_phase_space() { return particle_data_id % 2 ? particles1 : particles0; }
	inline DevicePlanetPhaseSpace& planet_phase_space() { return planet_data_id % 2 ? planets1 : planets0; }
};

struct HostParticlePhaseSpace
{
	size_t n, n_alive;

	Hvf64_3 r, v, rj, vj;
	Hvf64_3 a;

	Hvu8 flags;
	Hvf32 deathtime;
	Hvu32 deathtime_index;
	Hvu32 id;

	inline HostParticlePhaseSpace() { }
	inline HostParticlePhaseSpace(size_t n) : n(n), n_alive(n), r(n), v(n), a(n), flags(n), deathtime(n), deathtime_index(n), id(n) { }
};

struct HostPlanetPhaseSpace
{
	size_t n, n_alive;
	Hvf64 m, eta;
	Hvf64_3 r, v, rj, vj;
	Hvf64_3 a;

	Hvf64_3 r_log;
	Hvf64_3 h0_log;

	f64_3 bary_r, bary_v;

	Hvu8 flags;
	Hvf32 deathtime;
	Hvu32 deathtime_index;
	Hvu32 id;

	inline HostPlanetPhaseSpace() { }
	inline HostPlanetPhaseSpace(size_t n, size_t tb_size) :
		n(n), n_alive(n), m(n), eta(n), r(n), v(n), rj(n), vj(n),
		a(n),
		r_log(n * tb_size), h0_log(tb_size),
       		flags(n), deathtime(n), deathtime_index(n), id(n) { }
};

struct HostData
{
	HostParticlePhaseSpace particles, encounter_particles;
	HostPlanetPhaseSpace planets;

	size_t tbsize;

	inline HostData() { }
};


bool load_data(HostData& hd, std::string plin, std::string icsin, size_t tbsize, size_t max_particle = 0, bool readmomenta = true);
void save_data(const HostData& hd, std::string plout, std::string icsout);
