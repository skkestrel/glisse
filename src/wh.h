#pragma once
#include "data.h"
#include "util.h"

#include <unordered_map>

namespace sr
{
	namespace wh
	{
		using namespace sr::data;

		bool kepeq(double dM, double ecosEo, double esinEo, double *dE, double *sindE, double *cosdE, uint32_t *iterations);
		bool kepeq_fixed(double dM, double ecosEo, double esinEo, double *dE, double *sindE, double *cosdE, uint32_t iterations);

		class WHIntegrator
		{
		public:
			Vf64 planet_inverse_helio_cubed, planet_inverse_jacobi_cubed;

			Vf64 particle_dist, particle_energy, particle_vdotr;
			Vf64 planet_dist, planet_energy, planet_vdotr;

			Vu8 particle_mask;
			Vu8 planet_mask;

			Vf64 planet_mu, planet_eta;
			Vf64 particle_mu;

			Vf64_3 planet_rj, planet_vj;
			Vf64_3 planet_a, particle_a;

			sr::util::LogQuartet<Vf64_3> planet_h0_log;

			Vf64 planet_rh;
			double inner_bound, outer_bound, hill_factor;

			size_t tbsize;

			double dt;

			WHIntegrator();
			WHIntegrator(HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, const Configuration &config);

			void recalculate_rh(const HostPlanetPhaseSpace &pl);
			void swap_logs();

			void integrate_planets_timeblock(HostPlanetPhaseSpace &pl, float64_t t);
			void integrate_particles_timeblock(const HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, size_t begin, size_t length, float64_t t);
			void gather_particles(const std::vector<size_t> &indices, size_t begin, size_t length);

			void step_planets(HostPlanetPhaseSpace &pl, float64_t t, size_t timestep_index);
			void step_particles(const HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, size_t begin, size_t length, float64_t t, size_t timestep_index);

			static bool drift_single(float64_t t, float64_t mu, f64_3 *r, f64_3 *v);
			static void drift(float64_t t, Vf64_3 &r, Vf64_3 &v, size_t start, size_t n, Vf64 &dist, Vf64 &energy, Vf64 &vdotr, Vf64 &mu, Vu8 &mask);

			template <bool old>
			void helio_acc_particle(const HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, size_t particle_index, float64_t time, size_t timestep_index);

			template <bool old>
			void helio_acc_particles(const HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &p, size_t begin, size_t length, float64_t time, size_t timestep_index);

			void helio_acc_planets(HostPlanetPhaseSpace &p, size_t index);
		};

		void calculate_planet_metrics(const HostPlanetPhaseSpace &p, double *energy, f64_3 *l);
	} // namespace wh
} // namespace sr
