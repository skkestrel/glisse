#include "wh.h"
#include "convert.h"

#include <iomanip>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace sr
{
	namespace wh
	{
		const uint32_t MAXKEP = 10;
		const float64_t TOLKEP = 1E-14;

		using namespace sr::data;

		void print_tiss(const HostPlanetPhaseSpace &pl, const HostParticlePhaseSpace &pa)
		{
			double aout, eout, iout, aj;
			sr::convert::to_elements(pl.m()[0] + pl.m()[1], pl.r()[1], pl.v()[1], nullptr, &aj);

			for (uint32_t k = 0; k < pa.n(); k++)
			{
				sr::convert::to_elements(pl.m()[0], pa.r()[k], pa.v()[k], nullptr, &aout, &eout, &iout);
				std::cerr << pa.id()[k] << " " << aj / aout + 2 * std::sqrt((1 - eout * eout) * aout / aj) * std::cos(iout) << std::endl;
			}
		}

		bool kepeq(double dM, double esinEo, double ecosEo, double *dE, double *sindE, double *cosdE, uint32_t *iterations)
		{
			double f, fp, delta;

			*sindE = std::sin(*dE);
			*cosdE = std::cos(*dE);

			for (uint32_t i = 0; i < MAXKEP; i++)
			{
				f = *dE - ecosEo * (*sindE) + esinEo * (1. - *cosdE) - dM;
				fp = 1. - ecosEo * (*cosdE) + esinEo * (*sindE);
				delta = -f / fp;
				if (std::fabs(delta) < TOLKEP)
				{
					*iterations = i;
					return false;
				}

				*dE += delta;
				*sindE = std::sin(*dE);
				*cosdE = std::cos(*dE);
			}

			return true;
		}

		void calculate_planet_metrics(const HostPlanetPhaseSpace &pl, double *energy, f64_3 *l)
		{
			f64_3 bary_r, bary_v;
			double bary_m;

			sr::convert::find_barycenter(pl.r(), pl.v(), pl.m(), pl.n_alive(), bary_r, bary_v, bary_m);

			Vf64_3 r(pl.n_alive()), v(pl.n_alive());
			for (size_t i = 0; i < pl.n_alive(); i++)
			{
				r[i] = pl.r()[i] - bary_r;
				v[i] = pl.v()[i] - bary_v;
			}

			if (energy)
			{
				double ke = 0.0;
				double pe = 0.0;

				for (size_t i = 0; i < pl.n_alive(); i++)
				{
					ke += 0.5 * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z) * pl.m()[i];
				}

				for (size_t i = 0; i < pl.n_alive() - 1; i++)
				{
					for (size_t j = i + 1; j < pl.n_alive(); j++)
					{
						double dx = r[i].x - r[j].x;
						double dy = r[i].y - r[j].y;
						double dz = r[i].z - r[j].z;

						pe -= pl.m()[i] * pl.m()[j] / std::sqrt(dx * dx + dy * dy + dz * dz);
					}
				}

				*energy = ke + pe;
			}

			if (l)
			{
				*l = f64_3(0);

				for (size_t i = 0; i < pl.n_alive(); i++)
				{
					*l += r[i].cross(v[i]) * pl.m()[i];
				}
			}
		}

		WHIntegrator::WHIntegrator() {}
		WHIntegrator::WHIntegrator(HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, const Configuration &config)
		{
			planet_inverse_helio_cubed = planet_inverse_jacobi_cubed = Vf64(pl.n());
			planet_dist = planet_energy = planet_vdotr = Vf64(pl.n());
			particle_dist = particle_energy = particle_vdotr = Vf64(pa.n());

			planet_mu = Vf64(pl.n());
			particle_mu = Vf64(pa.n());

			planet_mask = Vu8(pl.n());
			particle_mask = Vu8(pa.n());

			planet_eta = Vf64(pl.n());
			planet_rj = planet_vj = Vf64_3(pl.n());
			planet_a = Vf64_3(pl.n());
			particle_a = Vf64_3(pa.n());
			planet_rh = Vf64(pl.n());

			tbsize = config.tbsize;
			dt = config.dt;

			outer_bound = config.outer_bound;
			hill_factor = config.hill_factor;

			planet_h0_log = sr::util::LogQuartet<Vf64_3>(tbsize);

			// Planet index 0 is always the sun.
			planet_eta[0] = pl.m()[0];

			for (size_t i = 1; i < pl.n(); i++)
			{
				planet_eta[i] = planet_eta[i - 1] + pl.m()[i];
				// a particle will be "killed" when it comes within planet_rh[i] of planet i. One might want to set this number
				// to be a multiple of the planet's respective hill sphere.
			}

			// [Probe] Set inner_bound as the radius of the Sun.
			planet_rh[0] = config.inner_bound;
			recalculate_rh(pl);

			sr::convert::helio_to_jacobi_r_planets(pl, planet_eta, planet_rj);
			sr::convert::helio_to_jacobi_v_planets(pl, planet_eta, planet_vj);

			std::copy(pl.r().begin() + 1, pl.r().end(), pl.r_log().old.begin());
			helio_acc_planets(pl, 0);
			helio_acc_particles<true>(pl, pa, 0, pa.n_alive(), 0, 0);
		}

		void WHIntegrator::recalculate_rh(const HostPlanetPhaseSpace &pl)
		{
			// for the rest of the planets, use the factor times the hill sphere
			for (size_t i = 1; i < pl.n(); i++)
			{
				double a, e;
				sr::convert::to_elements(pl.m()[0] + pl.m()[i], pl.r()[i], pl.v()[i], nullptr, &a, &e);
				planet_rh[i] = hill_factor * a * (1 - e) * std::pow(pl.m()[i] / (3 * pl.m()[0]), 1. / 3);
			}
		}

		void WHIntegrator::swap_logs()
		{
			planet_h0_log.swap_logs();
		}

		void WHIntegrator::integrate_planets_timeblock(HostPlanetPhaseSpace &pl, float64_t t)
		{
			for (size_t i = 0; i < tbsize; i++)
			{
				step_planets(pl, t, i);
				t += dt;
			}
		}

		void WHIntegrator::integrate_particles_timeblock(const HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, size_t begin, size_t length, float64_t t)
		{
			if (pa.deathtime_index().empty())
			{
				throw std::invalid_argument("Deathtime index array not allocated");
			}

			for (size_t i = begin; i < begin + length; i++)
			{
				this->particle_mu[i] = pl.m()[0];
			}

			for (size_t i = 0; i < tbsize; i++)
			{
				step_particles(pl, pa, begin, length, t, i);
				t += dt;
			}
		}

		template <bool old>
		void WHIntegrator::helio_acc_particle(const HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, size_t particle_index, float64_t time, size_t timestep_index)
		{
			f64_3 &a = particle_a[particle_index];
			a = planet_h0_log.get<old>()[timestep_index];

			for (size_t j = 1; j < pl.n_alive(); j++)
			{
				f64_3 dr = pa.r()[particle_index] - pl.r_log().get<old>()[pl.log_index_at<old>(timestep_index, j)];
#ifdef USE_FMA
				float64_t planet_rji2 = std::fma(dr.x, dr.x, std::fma(dr.y, dr.y, dr.z * dr.z));
#else
				float64_t planet_rji2 = dr.lensq();
#endif

				float64_t irij3 = 1. / (planet_rji2 * std::sqrt(planet_rji2));
				float64_t fac = pl.m()[j] * irij3;

#ifdef USE_FMA
				a.x = std::fma(-dr.x, fac, a.x);
				a.y = std::fma(-dr.y, fac, a.y);
				a.z = std::fma(-dr.z, fac, a.z);
#else
				a -= dr * fac;
#endif

				if (planet_rji2 < planet_rh[j] * planet_rh[j])
				{
					pa.deathflags()[particle_index] = pa.deathflags()[particle_index] & 0x00FF;
					pa.deathflags()[particle_index] = static_cast<uint16_t>(pa.deathflags()[particle_index] | (j << 8) | 0x0001);
				}
			}

			float64_t planet_rji2 = pa.r()[particle_index].lensq();

			if (planet_rji2 < planet_rh[0] * planet_rh[0])
			{
				pa.deathflags()[particle_index] = pa.deathflags()[particle_index] & 0x00FF;
				pa.deathflags()[particle_index] = pa.deathflags()[particle_index] | 0x0001;
			}

			if (planet_rji2 > outer_bound * outer_bound)
			{
				pa.deathtime()[particle_index] = static_cast<float>(time);
				pa.deathflags()[particle_index] = pa.deathflags()[particle_index] | 0x0002;
			}
		}

		template <bool old>
		void WHIntegrator::helio_acc_particles(const HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, size_t begin, size_t length, float64_t time, size_t timestep_index)
		{
			for (size_t i = begin; i < begin + length; i++)
			{
				helio_acc_particle<old>(pl, pa, i, time, timestep_index);
			}
		}

		void WHIntegrator::helio_acc_planets(HostPlanetPhaseSpace &p, size_t index)
		{
			for (size_t i = 1; i < p.n_alive(); i++)
			{
				float64_t r2 = p.r()[i].lensq();
				this->planet_inverse_helio_cubed[i] =
					1. / (std::sqrt(r2) * r2);
				r2 = this->planet_rj[i].lensq();
				this->planet_inverse_jacobi_cubed[i] =
					1. / (std::sqrt(r2) * r2);
			}

			// compute common heliocentric acceleration
			f64_3 a_common(0);
			for (size_t i = 2; i < p.n_alive(); i++)
			{
				float64_t mfac = p.m()[i] * this->planet_inverse_helio_cubed[i];
				a_common -= p.r()[i] * mfac;
			}

			// Load this into all the arrays
			for (size_t i = 1; i < p.n_alive(); i++)
			{
				planet_a[i] = a_common;
			}

			planet_h0_log.get<true>()[index] = a_common - p.r()[1] * p.m()[1] * this->planet_inverse_helio_cubed[1];

			// Now do indirect acceleration ; note that planet 1 does not receive a contribution
			for (size_t i = 2; i < p.n_alive(); i++)
			{
				planet_a[i] += (this->planet_rj[i] * this->planet_inverse_jacobi_cubed[i] - p.r()[i] * this->planet_inverse_helio_cubed[i]) * p.m()[0];
			}

			/* next term ; again, first planet does not participate */
			f64_3 a_accum(0);
			for (size_t i = 2; i < p.n_alive(); i++)
			{
				float64_t mfac = p.m()[i] * p.m()[0] * this->planet_inverse_jacobi_cubed[i] / this->planet_eta[i - 1];
				a_accum += this->planet_rj[i] * mfac;
				planet_a[i] += a_accum;
			}

			/* Finally, incorporate the direct accelerations */
			for (size_t i = 1; i < p.n_alive() - 1; i++)
			{
				for (size_t j = i + 1; j < p.n_alive(); j++)
				{
					f64_3 dr = p.r()[j] - p.r()[i];
					float64_t r2 = dr.lensq();
					float64_t irij3 = 1. / (r2 * std::sqrt(r2));

					float64_t mfac = p.m()[i] * irij3;
					planet_a[j] -= dr * mfac;

					// acc. on i is just negative, with m[j] instead
					mfac = p.m()[j] * irij3;
					planet_a[i] += dr * mfac;
				}
			}
		}

		bool WHIntegrator::drift_single(float64_t dt, float64_t mu, f64_3 *r, f64_3 *v)
		{
			float64_t dist, vsq, vdotr;
			dist = std::sqrt(r->lensq());
			vsq = v->lensq();
			vdotr = v->x * r->x + v->y * r->y + v->z * r->z;

			float64_t energy = vsq;
			energy *= 0.5;
			energy -= mu / dist;

			if (energy < 0)
			{
				float64_t a = -0.5 * mu / energy;
				float64_t n_ = std::sqrt(mu / (a * a * a));
				float64_t ecosEo = 1.0 - dist / a;
				float64_t esinEo = vdotr / (n_ * a * a);
				// Unused variable esq:
				// float64_t esq = ecosEo * ecosEo + esinEo * esinEo;

				// subtract off an integer multiple of complete orbits
				float64_t dM = dt * n_ - M_2PI * (int)(dt * n_ / M_2PI);

				// remaining time to advance
				dt = dM / n_;
				// maybe parallelize this

				// call kepler equation solver with initial guess in dE already
				float64_t dE = dM - esinEo + esinEo * std::cos(dM) + ecosEo * std::sin(dM);
				float64_t sindE, cosdE;

				uint32_t its;
				if (kepeq(dM, esinEo, ecosEo, &dE, &sindE, &cosdE, &its))
				{
					return true;
				}

				float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
				float64_t f = 1.0 + a * (cosdE - 1.0) / dist;
				float64_t g = dt + (sindE - dE) / n_;
				float64_t fdot = -n_ * sindE * a / (dist * fp);
				float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

				f64_3 r0 = *r;
				f64_3 v0 = *v;
				*r = r0 * f + v0 * g;
				*v = r0 * fdot + v0 * gdot;
				return false;
			}
			// hyperbolic orbit...
			return true;
		}

		void WHIntegrator::drift(float64_t t, Vf64_3 &r, Vf64_3 &v, size_t start, size_t n, Vf64 &dist, Vf64 &energy, Vf64 &vdotr, Vf64 &mu, Vu8 &mask)
		{
			for (size_t i = start; i < start + n; i++)
			{
				dist[i] = std::sqrt(r[i].lensq());
				energy[i] = v[i].lensq();
				vdotr[i] = v[i].x * r[i].x + v[i].y * r[i].y + v[i].z * r[i].z;
			}

			for (size_t i = start; i < start + n; i++)
			{
				energy[i] *= 0.5;
				energy[i] -= mu[i] / dist[i];
			}

			for (size_t i = start; i < start + n; i++)
			{
				if (mask[i])
					continue;
				if (energy[i] >= 0)
				{
					std::ostringstream ss;
					ss << "unbound orbit of planet " << i << " energy = " << energy[i] << std::endl;

					/*
				for (size_t j = start; j < start + n; j++)
				{
					ss << "p " << r[j].x << " " << r[j].y << " " << r[j].z << std::endl;
					ss << "v " << v[j].x << " " << v[j].y << " " << v[j].z << std::endl;
				}
				*/

					throw std::runtime_error(ss.str());
				}
				else
				{
					f64_3 r0 = r[i];
					f64_3 v0 = v[i];

					// maybe parallelize this
					float64_t a = -0.5 * mu[i] / energy[i];
					float64_t n_ = std::sqrt(mu[i] / (a * a * a));
					float64_t ecosEo = 1.0 - dist[i] / a;
					float64_t esinEo = vdotr[i] / (n_ * a * a);
					// float64_t e = std::sqrt(ecosEo * ecosEo + esinEo * esinEo);

					// subtract off an integer multiple of complete orbits
					float64_t dM = t * n_ - M_2PI * (int)(t * n_ / M_2PI);

					// remaining time to advance
					float64_t adv_dt = dM / n_;

					// call kepler equation solver with initial guess in dE already
					float64_t dE = dM - esinEo + esinEo * std::cos(dM) + ecosEo * std::sin(dM);
					float64_t sindE, cosdE;

					bool error;
					uint32_t its;
					error = kepeq(dM, esinEo, ecosEo, &dE, &sindE, &cosdE, &its);

					if (error)
					{
						throw std::runtime_error("Unconverging kepler");
					}

					float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
					float64_t f = 1.0 + a * (cosdE - 1.0) / dist[i];
					float64_t g = adv_dt + (sindE - dE) / n_;
					float64_t fdot = -n_ * sindE * a / (dist[i] * fp);
					float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

					r[i] = r0 * f + v0 * g;
					v[i] = r0 * fdot + v0 * gdot;
				}
			}
		}

		void WHIntegrator::step_particles(const HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, size_t begin, size_t length, float64_t t, size_t timestep_index)
		{
			for (size_t i = begin; i < begin + length; i++)
			{
				this->particle_mask[i] = pa.deathflags()[i] != 0;

				if (!this->particle_mask[i])
				{
					pa.v()[i] += this->particle_a[i] * (dt / 2);
				}
			}

			// Drift all the particles along their Jacobi Kepler ellipses
			// Can change false to true to use fixed iterations
			drift(dt, pa.r(), pa.v(), begin, length, particle_dist, particle_energy, particle_vdotr, particle_mu, particle_mask);

			// find the accelerations of the heliocentric velocities
			helio_acc_particles<false>(pl, pa, begin, length, t, timestep_index);

			for (size_t i = begin; i < begin + length; i++)
			{
				if (!this->particle_mask[i])
				{
					pa.v()[i] += particle_a[i] * (dt / 2);
					pa.deathtime_index()[i] = static_cast<uint32_t>(timestep_index + 1);
				}
			}
		}

		void WHIntegrator::gather_particles(const std::vector<size_t> &indices, size_t begin, size_t length)
		{
			gather(particle_a, indices, begin, length);
		}

		void WHIntegrator::step_planets(HostPlanetPhaseSpace &pl, float64_t t, size_t timestep_index)
		{
			// std::cerr << "pl. " << t << " " << pl.r()[1] << " " << pl.v()[1] << std::endl;

			(void)t;

			double new_dt = dt;

			for (size_t i = 1; i < pl.n_alive(); i++)
			{
				pl.v()[i] += planet_a[i] * (new_dt / 2);
			}

			// Convert the heliocentric velocities to Jacobi velocities
			sr::convert::helio_to_jacobi_v_planets(pl, planet_eta, planet_vj);

			for (size_t i = 1; i < pl.n_alive(); i++)
			{
				// Each Jacobi Kepler problem has a different central mass
				this->planet_mu[i] = pl.m()[0] * this->planet_eta[i] / this->planet_eta[i - 1];
				this->planet_mask[i] = 0;
			}

			// Drift all the particles along their Jacobi Kepler ellipses
			drift(new_dt, this->planet_rj, this->planet_vj, 1, pl.n_alive() - 1, planet_dist, planet_energy, planet_vdotr, planet_mu, planet_mask);

			// convert Jacobi vectors to helio. ones for acceleration calc
			sr::convert::jacobi_to_helio_planets(planet_eta, planet_rj, planet_vj, pl);

			// find the accelerations of the heliocentric velocities
			helio_acc_planets(pl, timestep_index);

			Vf64_3 *r_log = &pl.r_log().old;
			Vf64_3 *v_log = &pl.v_log().old;

			std::copy(pl.r().begin() + 1, pl.r().end(), r_log->begin() + (pl.n_alive() - 1) * timestep_index);
			std::copy(pl.v().begin() + 1, pl.v().end(), v_log->begin() + (pl.n_alive() - 1) * timestep_index);

			for (size_t i = 1; i < pl.n_alive(); i++)
			{
				pl.v()[i] += planet_a[i] * (new_dt / 2);
			}
		}
	} // namespace wh
} // namespace sr
