#include "convert.h"

// requires pl to have bary_r calculated
void helio_to_jacobi_r_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p)
{
	for (size_t i = 0; i < p.n; i++)
	{
		p.rj[i] = p.r[i] - pl.bary_r;
	}
}

void helio_to_jacobi_v_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p)
{
	for (size_t i = 0; i < p.n; i++)
	{
		p.vj[i] = p.v[i] - pl.bary_v;
	}
}

void helio_to_jacobi_v_planets(HostPlanetPhaseSpace& p)
{
	// COM
	p.vj[0] = f64_3(0);

	// same as heliocentric
	p.vj[1] = p.v[1];

	// momentum sum
	f64_3 psum(0);

	for (size_t i = 2; i < p.n; i++)
	{
		// velocity of interior COM
		f64_3 vsum;

		psum += p.v[i - 1] * p.m[i - 1];
		vsum = psum / p.eta[i - 1];

		p.vj[i] = p.v[i] - vsum;
	}

	psum += p.v[p.n - 1] * p.m[p.n - 1];
	p.bary_v = psum / p.eta[p.n - 1];
}

void helio_to_jacobi_r(HostPlanetPhaseSpace& p)
{
	p.eta[0] = p.m[0];

	for (size_t i = 1; i < p.n; i++)
	{
		p.eta[i] = p.eta[i - 1] + p.m[i];
	}

	// pick origin at baricenter
	p.rj[0] = f64_3(0);

	// first jacobi coordinate is same as heliocentric
	p.rj[1] = p.r[1];

	f64_3 sum(0);
	f64_3 bary;

	for (size_t i = 2; i < p.n; i++)
	{
		sum += p.r[i - 1] * p.m[i - 1];
		bary = sum / p.eta[i - 1];

		p.rj[i] = p.r[i] - bary;
	}

	sum += p.r[p.n - 1] * p.m[p.n - 1];
	p.bary_r = sum / p.eta[p.n - 1];
}

void helio_to_bary(HostData& hd)
{
	f64_3 r(0);
	f64_3 v(0);

	double totalm = 0;

	for (size_t i = 0; i < hd.planets.n; i++)
	{
		r += hd.planets.r[i] * hd.planets.m[i];
		v += hd.planets.v[i] * hd.planets.m[i];

		totalm += hd.planets.m[i];
	}

	r /= totalm;
	v /= totalm;

	for (size_t i = 0; i < hd.planets.n; i++)
	{
		hd.planets.r[i] -= r;
		hd.planets.v[i] -= v;
	}

	for (size_t i = 0; i < hd.particles.n; i++)
	{
		hd.particles.r[i] -= r;
		hd.particles.v[i] -= v;
	}
}

void bary_to_helio(HostData& hd)
{
	f64_3 r;
	f64_3 v;
	r = hd.planets.r[0];
	v = hd.planets.v[0];

	for (size_t i = 0; i < hd.planets.n; i++)
	{
		hd.planets.r[i] -= r;
		hd.planets.v[i] -= v;
	}

	hd.planets.bary_r = -r;
	hd.planets.bary_v = -v;

	for (size_t i = 0; i < hd.particles.n; i++)
	{
		hd.particles.r[i] -= r;
		hd.particles.v[i] -= v;
	}
}
