#include "data.h"
#include <fstream>
#include <iomanip>
#include <limits>

bool load_data(HostData& hd, std::string plin, std::string icsin, size_t tbsize, size_t max_particle, bool readmomenta)
{
	hd.tbsize = tbsize;
	std::ifstream plinfile(plin), icsinfile(icsin);

	size_t npl;
	plinfile >> npl;

	hd.planets = HostPlanetPhaseSpace(npl, tbsize);

	for (size_t i = 0; i < npl; i++)
	{
		plinfile >> hd.planets.m[i];
		plinfile >> hd.planets.r[i].x >> hd.planets.r[i].y >> hd.planets.r[i].z;
		plinfile >> hd.planets.v[i].x >> hd.planets.v[i].y >> hd.planets.v[i].z;

		if (readmomenta)
		{
			hd.planets.v[i].x /= hd.planets.m[i];
			hd.planets.v[i].y /= hd.planets.m[i];
			hd.planets.v[i].z /= hd.planets.m[i];
		}
	}

	size_t npart;
	icsinfile >> npart;
	if (max_particle > 0) npart = std::min(npart, max_particle);

	hd.particles = HostParticlePhaseSpace(npart);

	for (size_t i = 0; i < npart; i++)
	{
		icsinfile >> hd.particles.r[i].x >> hd.particles.r[i].y >> hd.particles.r[i].z;
		icsinfile >> hd.particles.v[i].x >> hd.particles.v[i].y >> hd.particles.v[i].z;

		std::string s;
		icsinfile >> s;
		if (!isdigit(s[0]))
		{
			icsinfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			hd.particles.deathtime[i] = 0;
			hd.particles.id[i] = i;
			hd.particles.flags[i] = 0;
		}
		else
		{
			hd.particles.deathtime[i] = std::stod(s);
			icsinfile >> hd.particles.flags[i] >> hd.particles.id[i];
		}
	}

	return false;
}

void save_data(const HostData& hd, std::string plout, std::string icsout)
{
	std::ofstream ploutfile(plout), icsoutfile(icsout);

	ploutfile << hd.planets.n << std::endl;
	ploutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.planets.n; i++)
	{
		ploutfile << hd.planets.m[i] << std::endl;
		ploutfile << hd.planets.r[i].x << " " << hd.planets.r[i].y << " " << hd.planets.r[i].z << std::endl;
		ploutfile << hd.planets.v[i].x << " " << hd.planets.v[i].y << " " << hd.planets.v[i].z << std::endl;
	}

	icsoutfile << hd.particles.n << std::endl;
	icsoutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.particles.n; i++)
	{
		icsoutfile << hd.particles.r[i].x << " " << hd.particles.r[i].y << " " << hd.particles.r[i].z << std::endl;
		icsoutfile << hd.particles.v[i].x << " " << hd.particles.v[i].y << " " << hd.particles.v[i].z << std::endl;
		icsoutfile << hd.particles.deathtime[i] << " " << hd.particles.flags[i] << " " << hd.particles.id[i] << std::endl;
	}
}
