#include "../src/data.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

static const char USAGE[] = R"(track-info
Usage:
    track-info [options] <input>

Options:
    -h, --help                     Show this screen.
)";

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "track-info");

	try
	{
		std::string inpath = args["<input>"].asString();

		sr::data::TrackReaderOptions opt;
		opt.take_all_particles = true;
		opt.take_all_planets = true;
		opt.silent = true;

		double lasttime;
		int n = 0;

		sr::data::read_tracks(inpath, opt,
			[&](sr::data::HostPlanetSnapshot& pl, sr::data::HostParticleSnapshot& pa, double time)
			{
				if (n == 0)
				{
					n++;
					std::cout << "There are " << pl.n << " planets:" << std::endl;
					for (size_t i = 0; i < pl.n; i++)
					{
						std::cout << pl.id[i] << ", ";
					}
					std::cout << std::endl;

					std::cout << "There are " << pa.n << " particles:" << std::endl;
					for (size_t i = 0; i < std::min(static_cast<size_t>(20), pa.n); i++)
					{
						std::cout << pa.id[i] << ", ";
					}
					if (pa.n > 20)
					{
						std::cout << "...";
					}
					std::cout << std::endl;
				}
				if (n == 1)
				{
					n++;
					std::cout << "The time step is " << time - lasttime << std::endl;
				}

				opt.take_all_particles = false;
				opt.take_all_planets = false;
				lasttime = time;
			});
		std::cout << "The length is " << lasttime << std::endl;
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
