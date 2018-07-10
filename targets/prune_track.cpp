#include "../data.h"
#include "../util.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include "../cxxopts.h"
#pragma GCC diagnostic pop

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

int main(int argc, char** argv)
{
	cxxopts::Options options("prune-track", "Prune large track files generated by simulations");
	options.add_options()
		("i,input", "Input directory", cxxopts::value<std::string>())
		("o,output", "Output directory", cxxopts::value<std::string>())
		("s,split-output", "Split output every n bytes", cxxopts::value<uint64_t>())
		("P,no-planets", "Throw away planets")
		("w,watch-particles", "Watch particles (comma separated list)", cxxopts::value<std::string>());

	auto result = options.parse(argc, argv);

	if (result.count("i") == 0)
	{
		std::cout << "Required flag -i" << std::endl;
		return -1;
	}
	if (result.count("o") == 0)
	{
		std::cout << "Required flag -i" << std::endl;
		return -1;
	}

	std::vector<uint32_t> particles;
	if (result.count("w") > 0)
	{
		std::istringstream ss(result["w"].as<std::string>());
		std::string token;

		while (std::getline(ss, token, ','))
		{
			particles.push_back(static_cast<uint32_t>(std::stoul(token)));
		}
	}

	std::sort(particles.begin(), particles.end());

	uint64_t splitbytes = 0;
	if (result.count("s") > 0)
	{
		splitbytes = result["s"].as<uint64_t>();
	}

	bool killplanets = false;
	if (result.count("P") > 0)
	{
		killplanets = true;
	}

	std::string inpath = result["i"].as<std::string>();
	std::string outpath = result["o"].as<std::string>();

	std::ostringstream ss;

	ss << outpath;
	if (outpath[outpath.size() - 1] != '/') ss << '/';
	ss << "track.0.out";

	size_t outnum = 1;
	std::ofstream outfile(ss.str(), std::ios_base::binary);

	for (size_t i = 0; true; i++)
	{
		ss = std::ostringstream();
		ss << inpath;
		if (inpath[inpath.size() - 1] != '/') ss << '/';
		ss << "track." << i << ".out";

		if (!sr::util::does_file_exist(ss.str()))
		{
			if (i == 0)
			{
				std::cout << "Error: empty directory" << std::endl;
				return -1;
			}
			else
			{
				std::cout << i << " files read" << std::endl;
				return 0;
			}
		}

		std::ifstream input(ss.str(), std::ios_base::binary);

		while (input)
		{
			double time;
			sr::data::HostParticleSnapshot pa;
			sr::data::HostPlanetSnapshot pl;
			sr::data::load_elements(input, pl, pa, time);

			if (killplanets)
			{
				pl.n_alive = pl.n = 0;
			}

			sr::data::HostParticleSnapshot paout;
			size_t lowbound = 0;
			for (size_t j = 0; i < particles.size(); i++)
			{
				size_t index = std::lower_bound(pa.id.begin() + lowbound, pa.id.end(), particles[j]) - pa.id.begin();
				lowbound = index;

				if (index != pa.n && pa.id[index] == particles[j])
				{
					paout.r.push_back(pa.r[index]);
					paout.v.push_back(pa.v[index]);
					paout.id.push_back(pa.id[index]);
				}
			}
			paout.n_alive = paout.n = paout.r.size();

			sr::data::save_elements(outfile, pl, paout, time);

			if (outfile.tellp() > static_cast<int>(splitbytes))
			{
				ss = std::ostringstream();

				ss << outpath;
				if (outpath[outpath.size() - 1] != '/') ss << '/';
				ss << "track." << outnum++ << ".out";

				outfile = std::ofstream(ss.str(), std::ios_base::binary);
			}
		}
	}

	return 0;
}
