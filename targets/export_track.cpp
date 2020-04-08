#include "../src/data.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

static const char USAGE[] = R"(export-track
Usage:
    export-track [options] <input> <output>

Options:
    -h, --help                     Show this screen.
    --true-anomaly                 Export the last number on each line as the true anomaly instead of the mean anomaly.
    --precision <val>              Export with val digits of precision. [default: 5]
    --radian                       Export with radians instead of degrees.
)";

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "export-track");

	try
	{
		std::string inpath = args["<input>"].asString();
		std::string outpath = args["<output>"].asString();
		int precision = std::stoi(args["--precision"].asString());
		bool radian = static_cast<bool>(args["--radian"]);
		bool trueanomaly = static_cast<bool>(args["--true-anomaly"]);

		std::ofstream out(outpath, std::ios_base::binary);

		sr::data::TrackReaderOptions opt;
		opt.take_all_particles = true;
		opt.take_all_planets = true;

		double mul = 180 / M_PI;
		if (!radian) mul = 1;

		sr::data::read_tracks(inpath, opt,
			[&](sr::data::HostPlanetSnapshot& pl, sr::data::HostParticleSnapshot& pa, double time)
			{
				for (size_t i = 0; i < pl.n; i++)
				{
					out << std::setprecision(15);
					out << -static_cast<int>(pl.id[i]) << " " << time << " ";
					out << std::setprecision(precision);
					out << pl.r[i].x << " " << pl.r[i].y << " " << mul*pl.r[i].z << " ";
					out << mul*pl.v[i].x << " " << mul*pl.v[i].y << " ";
					if (trueanomaly)
					{
						double cosf = std::cos(pl.v[i].z);
						double E = std::acos((pl.r[i].y + cosf) / (1 + pl.r[i].y * cosf));
						// E = E * sign(f)
						E = pl.v[i].z > 0 ? E : -E;
						double M = E - pl.r[i].y * std::sin(E);
						out << mul*M << std::endl;
					}
					else out << mul*pl.v[i].z << std::endl;
					
				}
				for (size_t i = 0; i < pa.n; i++)
				{
					out << std::setprecision(15);
					out << pa.id[i] << " " << time << " ";
					out << std::setprecision(precision);
					out << pa.r[i].x << " " << pa.r[i].y << " " << mul*pa.r[i].z << " ";
					out << mul*pa.v[i].x << " " << mul*pa.v[i].y << " ";
					if (trueanomaly)
					{
						double cosf = std::cos(pa.v[i].z);
						double E = std::acos((pa.r[i].y + cosf) / (1 + pa.r[i].y * cosf));
						// E = E * sign(f)
						E = pa.v[i].z > 0 ? E : -E;
						double M = E - pa.r[i].y * std::sin(E);
						out << mul*M << std::endl;
					}
					else out << mul*pa.v[i].z << std::endl;
				}
			});
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
