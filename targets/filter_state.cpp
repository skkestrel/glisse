#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <iomanip>

#include "../docopt/docopt.h"
#include "../src/data.h"
#include "../src/wh.h"
#include "../src/convert.h"
#include "../src/util.h"

static const char USAGE[] = R"(filter-state
Usage:
    filter-state [options] <input> [<criteria>...]

Options:
    -h, --help                         Show this screen.
    -i <path>, --initial-state <path>  Initial state input file
    -c <list>, --csv <list>            Comma-separated list of paths to csvs that can contain attributes to filter by
    -u, --union                        Take union of criteria instead of intersection
    -b, --binary                       Read binary input
    -B, --barycentric                  Calculate barycentric elements
    -o <file>, --output <file>         Write output filtered state
)";

struct Criterion
{
	std::string variable;
	int comparison; // -1: var lt konst, 0: var eq konst, 1: var gt const
	double konst;
};

struct CsvFile
{
	std::vector<std::vector<double>> values;
	std::unordered_map<std::string, size_t> name_to_vector_index;
	std::unordered_map<uint32_t, size_t> id_to_index;
};

const double EPS = 1e-13;
int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "filter-state");

	try
	{
		bool use_union = args["--union"].asBool();

		std::vector<Criterion> criteria;

		if (args["<criteria>"])
		{
			for (const std::string& str : args["<criteria>"].asStringList())
			{
				Criterion crit;
				size_t foundgt = str.find('>');
				size_t foundeq = str.find('=');
				size_t foundlt = str.find('<');

				size_t split;

				if (foundgt != std::string::npos)
				{
					crit.comparison = 1;
					split = foundgt;
				}
				else if (foundeq != std::string::npos)
				{
					crit.comparison = 0;
					split = foundeq;
				}
				else if (foundlt != std::string::npos)
				{
					crit.comparison = -1;
					split = foundlt;
				}
				else
				{
					std::ostringstream ss;
					ss << "Could not parse criterion " << str;
					throw std::runtime_error(ss.str());
				}

				crit.variable = str.substr(0, split);
				crit.konst = std::stod(str.substr(split + 1, std::string::npos));

				criteria.push_back(crit);
			}	
		}

		std::vector<CsvFile> files;
		std::unordered_map<std::string, size_t> name_to_file_index;

		if (args["--csv"])
		{
			std::istringstream ss(args["--csv"].asString());
			std::string token;

			while (std::getline(ss, token, ','))
			{
				std::string filename = token;
				
				std::ifstream infile(filename);
				std::string line;
				size_t linenum = 0;
				CsvFile newfile;

				while (std::getline(infile, line))
				{
					std::istringstream ss2(line);
					std::string token2;
					size_t token_num = 0;

					while (std::getline(ss2, token2, ','))
					{
						if (linenum == 0)
						{
							newfile.name_to_vector_index[token2] = token_num;
							name_to_file_index[token2] = files.size();
						}
						else
						{
							if (newfile.name_to_vector_index["id"] == token_num)
							{
								newfile.id_to_index[static_cast<uint32_t>(std::stoull(token2))] = linenum - 1;
							}
							else
							{
								newfile.values[token_num].push_back(std::stod(token2));
							}
						}

						token_num++;
					}

					if (linenum == 0)
					{
						for (size_t i = 0; i < newfile.name_to_vector_index.size(); i++)
						{
							newfile.values.push_back(std::vector<double>());
						}
					}

					linenum++;
				}

				files.push_back(std::move(newfile));
			}
		}

		name_to_file_index.erase("id");

		sr::data::Configuration config = sr::data::Configuration::create_dummy();
		sr::data::HostData hd;
		sr::data::HostData hd_init;

		config.hybridin = args["<input>"].asString();
		config.readbinary = args["--binary"].asBool();

		load_data(hd.planets, hd.particles, config);
		hd.particles.sort_by_id(0, hd.particles.n());
		if (args["--barycentric"])
		{
			sr::convert::to_bary(hd);
		}
		else
		{
			sr::convert::to_helio(hd);
		}


		bool has_init = static_cast<bool>(args["--initial-state"]);
		if (has_init)
		{
			config.hybridin = args["--initial-state"].asString();
			load_data(hd_init.planets, hd_init.particles, config);

			hd_init.particles.sort_by_id(0, hd_init.particles.n());

			if (hd_init.particles.id() != hd.particles.id())
			{
				throw std::runtime_error("Initial state is not congruent with input state");
			}

			if (args["--barycentric"])
			{
				sr::convert::to_bary(hd_init);
			}
			else
			{
				sr::convert::to_helio(hd_init);
			}
		}

		std::vector<size_t> candidates;
		for (size_t i = 0; i < hd.particles.n(); i++)
		{
			int esign, esign_I;
			double A, E, I, CAPOM, OM, F;
			double A_I, E_I, I_I, CAPOM_I, OM_I, F_I;

			sr::convert::to_elements(hd.planets.m()[0], hd.particles.r()[i], hd.particles.v()[i], &esign, &A, &E, &I, &CAPOM, &OM, &F);

			if (has_init)
			{
				sr::convert::to_elements(hd_init.planets.m()[0], hd_init.particles.r()[i], hd_init.particles.v()[i], &esign_I, &A_I, &E_I, &I_I, &CAPOM_I, &OM_I, &F_I);
			}

			bool ok = !use_union;

			for (const auto& crit : criteria)
			{
				double val;
				if (crit.variable == "a") val = A;
				else if (crit.variable == "e") val = E;
				else if (crit.variable == "i") val = I;
				else if (crit.variable == "O") val = CAPOM;
				else if (crit.variable == "o") val = OM;
				else if (crit.variable == "f") val = F;
				else if (crit.variable == "id") val = static_cast<double>(hd.particles.id()[i]);
				else if (crit.variable == "a_i") val = A_I;
				else if (crit.variable == "e_i") val = E_I;
				else if (crit.variable == "i_i") val = I_I;
				else if (crit.variable == "O_i") val = CAPOM_I;
				else if (crit.variable == "o_i") val = OM_I;
				else if (crit.variable == "f_i") val = F_I;
				else if (crit.variable == "deathtime") val = hd.particles.deathtime()[i];
				else if (crit.variable == "killer")
				{
					if (hd.particles.deathflags()[i] == 0) val = -1;
					else val = static_cast<double>(hd.particles.deathflags()[i] >> 8);
				}
				else if (name_to_file_index.count(crit.variable) > 0)
				{
					const auto& file = files[name_to_file_index.at(crit.variable)];
					const auto& values = file.values[file.name_to_vector_index.at(crit.variable)];

					if (file.id_to_index.count(hd.particles.id()[i]) > 0)
					{
						val = values[file.id_to_index.at(hd.particles.id()[i])];
					}
					else
					{
						val = std::numeric_limits<double>::quiet_NaN();
					}
				}
				else
				{
					std::ostringstream ss;
					ss << "Unknown value " << crit.variable;
					throw std::runtime_error(ss.str());
				}

				bool newbool;
				if (crit.comparison == 1)
				{
					newbool = (val > crit.konst);
				}
				else if (crit.comparison == -1)
				{
					newbool = (val < crit.konst);
				}
				else if (crit.comparison == 0)
				{
					newbool = std::abs(val - crit.konst) < EPS;
				}
				else
				{
					throw std::runtime_error("Internal error");
				}

				if (use_union)
				{
					ok = ok || newbool;
					if (ok) break;
				}
				else
				{
					ok = ok && newbool;
					if (!ok) break;
				}
			}

			if (ok)
			{
				candidates.push_back(i);
			}
		}

		
		if (args["--output"])
		{
			sr::data::HostData hd_out;
			hd_out.planets = hd.planets;

			std::string outfile = args["--output"].asString();

			hd.particles.filter(candidates, hd_out.particles);

			sr::data::Configuration outcfg;
			outcfg.hybridout = outfile;
			sr::data::save_data(hd_out.planets.base, hd_out.particles, outcfg, outcfg.hybridout);
		}

		std::cout << "ID    | a        e        i       Om       om       f       ";
		for (auto& pair : name_to_file_index)
		{
			std::cout << std::left << std::setw(10) << pair.first;
		}
		std::cout << std::right << std::endl;

		int num = 0;
		for (auto& i : candidates)
		{
			int esign;
			double A, E, I, CAPOM, OM, F;
			sr::convert::to_elements(hd.planets.m()[0], hd.particles.r()[i], hd.particles.v()[i], &esign, &A, &E, &I, &CAPOM, &OM, &F);

			std::cout << std::setw(5) << hd.particles.id()[i] << " |";
			std::cout << std::fixed << std::setprecision(4)
				<< std::setw(8) << A << std::setw(8) << E << std::setw(8) << I
				<< std::setw(9) << CAPOM << std::setw(9) << OM << std::setw(9) << F << " ";
			for (auto& pair : name_to_file_index)
			{
				const auto& file = files[pair.second];
				if (file.id_to_index.count(hd.particles.id()[i]) > 0)
				{
					double val = file.values[file.name_to_vector_index.at(pair.first)][file.id_to_index.at(hd.particles.id()[i])];
					std::cout << std::setw(10) << std::defaultfloat << val;
				}
				else
				{
					std::cout << std::setw(10) << std::defaultfloat << std::numeric_limits<double>::quiet_NaN();
				}
			}
			std::cout << std::endl;
			num++;

			if (num % 10 == 0 && num != 0)
			{
				std::cout << candidates.size() - num << " left (" << candidates.size() << " total) (press return or q) ";
				int ch = '\0';
				while (ch != '\n' && ch != 'q')
				{
					ch = std::getchar();
				}

				if (ch == 'q') return 0;
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
