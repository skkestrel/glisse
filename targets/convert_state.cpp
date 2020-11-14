#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <iomanip>

#include "../src/data.h"
#include "../src/wh.h"
#include "../src/convert.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

static const char USAGE[] = R"(convert-state
Usage:
    convert-state [options] <commands>...

Options:
    -h, --help                     Show this screen.
)";

using namespace sr::data;
using namespace sr::convert;
const double EPS = 1e-13;

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "convert-state");

	try
	{
		Configuration config;
		HostData hd;

		bool ishelio;
		bool binary = false;
		bool momentum = false;

		const std::vector<std::string>& commands = args["<commands>"].asStringList();

		for (size_t i = 0; i < commands.size(); i++)
		{
			const std::string& arg = commands[i];

			if (arg == "momentum")
			{
				momentum = true;
			}
			else if (arg == "velocity")
			{
				momentum = false;
			}
			else if (arg == "binary")
			{
				binary = true;
			}
			else if (arg == "ascii")
			{
				binary = false;
			}
			else if (arg == "read")
			{
				config.hybridin = commands[i + 1];
				config.readbinary = binary;

				load_data(hd.planets, hd.particles, config);

				if (hd.planets.r()[0].lensq() < EPS && hd.planets.r()[0].lensq() < EPS) ishelio = true;
				else
				{
					to_bary(hd);
				}

				i++;
			}
			else if (arg == "read-split")
			{
				config.plin = commands[i + 1];
				config.icsin = commands[i + 2];

				config.readbinary = false;
				config.readsplit = true;

				load_data(hd.planets, hd.particles, config);

				if (hd.planets.r()[0].lensq() < EPS && hd.planets.r()[0].lensq() < EPS) ishelio = true;
				else
				{
					to_bary(hd);
				}

				i += 2;
			}
			else if (arg == "write-swift")
			{
				std::string plout = commands[i + 1];
				std::string icsout = commands[i + 2];
				std::ofstream o1(plout);
				std::ofstream o2(icsout);


				save_data_swift(hd.planets.base, hd.particles, o1, o2);

				i += 2;
			}
			else if (arg == "write")
			{
				config.hybridout = commands[i + 1];
				config.writebinary = binary;

				save_data(hd.planets.base, hd.particles, config, config.hybridout);

				i++;
			}
			else if (arg == "to-bary")
			{
				to_bary(hd);
			}
			else if (arg == "to-helio")
			{
				to_helio(hd);
			}
			else if (arg == "to-elements")
			{
				double a, e, I, capom, om, f;
				int esign;

				double totalmass = 0;
				for (size_t j = 0; j < hd.planets.n(); j++)
				{
					totalmass += hd.planets.m()[j];
				}

				for (size_t j = 1; j < hd.planets.n(); j++)
				{
					if (ishelio)
					{
						to_elements(hd.planets.m()[j] + hd.planets.m()[0], hd.planets.r()[j], hd.planets.v()[j], &esign, &a, &e, &I, &capom, &om, &f);
					}
					else
					{
						to_elements(totalmass, hd.planets.r()[j], hd.planets.v()[j], &esign, &a, &e, &I, &capom, &om, &f);
					}

					if (esign == 0)
					{
						std::cout << "Parabolic orbit detected!" << std::endl;
					}

					hd.planets.r()[j].x = a;
					hd.planets.r()[j].y = e;
					hd.planets.r()[j].z = I;
					hd.planets.v()[j].x = capom;
					hd.planets.v()[j].y = om;
					hd.planets.v()[j].z = f;
				}

				for (size_t j = 0; j < hd.particles.n(); j++)
				{
					if (ishelio)
					{
						to_elements(hd.planets.m()[0], hd.particles.r()[j], hd.particles.v()[j], &esign, &a, &e, &I, &capom, &om, &f);
					}
					else
					{
						to_elements(totalmass, hd.particles.r()[j], hd.particles.v()[j], &esign, &a, &e, &I, &capom, &om, &f);
					}

					if (esign == 0)
					{
						std::cout << "Parabolic orbit detected!" << std::endl;
					}

					hd.particles.r()[j].x = a;
					hd.particles.r()[j].y = e;
					hd.particles.r()[j].z = I;
					hd.particles.v()[j].x = capom;
					hd.particles.v()[j].y = om;
					hd.particles.v()[j].z = f;
				}
			}
			else
			{
				std::ostringstream ss;
				ss << "Unknown command " << arg;
				throw std::runtime_error(ss.str());
			}
		}
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
