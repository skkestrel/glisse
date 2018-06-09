#include "wh.h"
#include "convert.h"

#include <iostream>
#include <cmath>

const size_t MAXKEP = 10;
const float64_t TOLKEP = 1E-13;

void helio_acc_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p, float64_t time, size_t index)
{
	for (size_t i = 0; i < p.n; i++)
	{
		p.a[i] = pl.h0_log[index];
		for (size_t j = 1; j < pl.n; j++)
		{
			f64_3 dr = p.r[i] - pl.r[j];
			float64_t rji2 = dr.lensq();
			float64_t irij3 = 1. / (rji2 * std::sqrt(rji2));
			float64_t fac = pl.m[j] * irij3;

			p.a[i] -= dr * fac;

			if (rji2 < 0.5 * 0.5)
			{
				p.deathtime[i] = time;
				p.flags[i] = j;
			}
		}
	}
}

void helio_acc_planets(HostPlanetPhaseSpace& p, size_t index)
{
	Hvf64 inverse_helio_cubed(p.n), inverse_jacobi_cubed(p.n);

	for (size_t i = 1; i < p.n; i++)
	{
		float64_t r2 = p.r[i].lensq();
		inverse_helio_cubed[i] = 1. / (std::sqrt(r2) * r2);
		r2 = p.rj[i].lensq();
		inverse_jacobi_cubed[i] = 1. / (std::sqrt(r2) * r2);
        }
	
        // compute common heliocentric acceleration
	f64_3 a_common(0);
	for (size_t i = 2; i < p.n; i++)    
	{
		float64_t mfac = p.m[i] * inverse_helio_cubed[i];
		a_common -= p.r[i] * mfac;
        }

        // Load this into all the arrays
	for (size_t i = 1; i < p.n; i++)    
	{
		p.a[i] = a_common;
        }

	p.h0_log[index] = a_common - p.m[1] * inverse_helio_cubed[1];
	
	// Now do indirect acceleration ; note that planet 1 does not receive a contribution 
	for (size_t i = 2; i < p.n; i++)    
	{
		p.a[i] += (p.rj[i] * inverse_jacobi_cubed[i] - p.r[i] * inverse_helio_cubed[i]) * p.m[0];
        }
	
	/* next term ; again, first planet does not participate */
	f64_3 a_accum(0);
	for (size_t i = 2; i < p.n; i++)    
	{
		float64_t mfac = p.m[i] * p.m[0] * inverse_jacobi_cubed[i] / p.eta[i-1];
		a_accum += p.r[i] * mfac;
		p.a[i] += a_accum;
        }

	/* Finally, incorporate the direct accelerations */
	for (size_t i = 1; i < p.n - 1; i++)    
	{
		for (size_t j = i + 1; j < p.n; j++)    
		{
			f64_3 dr = p.r[j] - p.r[i];
			float64_t r2 = dr.lensq();
			float64_t irij3 = 1. / (r2 * std::sqrt(r2));

			float64_t mfac = p.m[i] * irij3;
			p.a[j] -= dr * mfac;

			// acc. on i is just negative, with m[j] instead
			mfac = p.m[j] * irij3;
			p.a[i] += dr * mfac;
		}
	}
}

size_t kepeq(double dM, double ecosEo, double esinEo, double* dE, double* sindE, double* cosdE)
{
	double f,fp, delta;

	*sindE = sin( *dE);
	*cosdE = cos( *dE);

	size_t i;
	for (i = 0; i < MAXKEP; i++)
	{
		f = *dE - ecosEo * (*sindE) + esinEo * (1. - *cosdE) - dM;
		fp = 1. - ecosEo * (*cosdE) + esinEo * (*sindE);
		delta = -f / fp;
		if (std::fabs(delta) < TOLKEP)
		{
			goto done;
		}

		*dE += delta;
		*sindE = std::sin(*dE);
		*cosdE = std::cos(*dE);
	}
	throw std::exception();

done:
	return i;
}

void drift(float64_t t, Hvu8& mask, Hvf64& mu, Hvf64_3& r, Hvf64_3& v, size_t start, size_t n)
{
	// save initial values
	Hvf64_3 r0 = r;
	Hvf64_3 v0 = v;

	Hvf64 dist(n);
	Hvf64 vsq(n);
	Hvf64 vdotr(n);
	for (size_t i = start; i < start + n; i++)
	{
		dist[i] = std::sqrt(r[i].lensq());
		vsq[i] = v[i].lensq();
		vdotr[i] = v0[i].x * r0[i].x + v0[i].y * r0[i].y + v0[i].z * r0[i].z;
	}

	Hvf64 energy = std::move(vsq);
	// vsq dies!
	for (size_t i = start; i < start + n; i++)
	{
		energy[i] *= 0.5;
		energy[i] -= mu[i] / dist[i];
	}

	for (size_t i = start; i < start + n; i++)
	{
		if (mask[i]) continue;
		if (energy[i] >= 0)
		{
			std::cerr << "unbound orbit" << std::endl;
			throw std::exception();
		}
		else
		{
			// maybe parallelize this
			float64_t a = -0.5 * mu[i] / energy[i];
			float64_t n = std::sqrt(mu[i] / (a * a * a));
			float64_t ecosEo = 1.0 - dist[i] / a;
			float64_t esinEo = vdotr[i] / (n * a * a);
			float64_t e = std::sqrt(ecosEo * ecosEo + esinEo * esinEo);

			// subtract off an integer multiple of complete orbits
			float64_t dM = t * n - M_2PI * (int) (t * n / M_2PI);

			// remaining time to advance
			float64_t dt = dM / n;

			// call kepler equation solver with initial guess in dE already
			float64_t dE = dM - esinEo + esinEo * std::cos(dM) + ecosEo * std::sin(dM);
			float64_t sindE, cosdE;
			kepeq(dM, ecosEo, esinEo, &dE, &sindE, &cosdE);

			float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
			float64_t f = 1.0 + a * (cosdE - 1.0) / dist[i];
			float64_t g = dt + (sindE - dE) / n;
			float64_t fdot = -n * sindE * a / (dist[i] * fp);
			float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

			r[i] = r0[i] * f + v0[i] * g;
			v[i] = r0[i] * fdot + v0[i] * gdot;
		}
	}
}

void initialize(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa)
{
	helio_to_jacobi_r_planets(pl);
	helio_to_jacobi_v_planets(pl);
}

void first_step(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, float64_t dt)
{
	helio_acc_planets(pl, 0);
	helio_acc_particles(pl, pa, 0, 0);

	for (size_t i = 1; i < pl.n; i++)
	{
		pl.v[i] += pl.a[i] * (dt / 2);
	}

	for (size_t i = 0; i < pa.n; i++)
	{
		pa.v[i] += pa.a[i] * (dt / 2);
	}
}

void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, float64_t t, size_t index, float64_t dt)
{
	Hvf64 mu(pa.n);
	Hvu8 mask(pa.n);
	for (size_t i = 0; i < pa.n; i++)
	{
		mu[i] = pl.m[0];
		mask[i] = pa.flags[i] > 0;
        }

	// Drift all the particles along their Jacobi Kepler ellipses
	drift(dt, mask, mu, pa.r, pa.v, 0, pa.n);

	// find the accelerations of the heliocentric velocities
	helio_acc_particles(pl, pa, t, index);

	for (size_t i = 0; i < pa.n; i++)
	{
		pa.v[i] += pa.a[i] * dt;
	}
}

void step_planets(HostPlanetPhaseSpace& pl, float64_t t, size_t index, float64_t dt)
{
	// Convert the heliocentric velocities to Jacobi velocities 
	helio_to_jacobi_v_planets(pl);

	Hvf64 mu(pl.n);
	Hvu8 mask(pl.n);
	for (size_t i = 1; i < pl.n; i++)
	{
		// Each Jacobi Kepler problem has a different central mass
		mu[i] = pl.m[0] * pl.eta[i] / pl.eta[i - 1];
		mask[i] = 0;
        }

	// Drift all the particles along their Jacobi Kepler ellipses
	drift(dt, mask, mu, pl.rj, pl.vj, 1, pl.n - 1);

	// convert Jacobi vectors to helio. ones for acceleration calc 
	jacobi_to_helio_planets(pl);

	// find the accelerations of the heliocentric velocities
	helio_acc_planets(pl, index);
	std::copy(pl.r.begin(), pl.r.end(), pl.r_log.begin() + (pl.n - 1) * index);

	for (size_t i = 1; i < pl.n; i++)
	{
		pl.v[i] += pl.a[i] * dt;
	}
}
