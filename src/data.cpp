#include "data.h"
#include "util.h"
#include "convert.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>

namespace sr
{
	namespace data
	{
		size_t stable_partition_alive_indices(const std::vector<uint16_t> &flags, size_t begin, size_t length, std::unique_ptr<std::vector<size_t>> *indices)
		{
			auto new_indices = std::make_unique<std::vector<size_t>>(length);
			std::iota(new_indices->begin(), new_indices->end(), 0);

			auto val = std::stable_partition(new_indices->begin(), new_indices->end(), [begin, &flags](size_t index) { return (flags[index + begin] & 0x00FE) == 0; }) - new_indices->begin() + begin;

			*indices = std::move(new_indices);
			return val;
		}

		size_t stable_partition_unflagged_indices(const std::vector<uint16_t> &flags, size_t begin, size_t length, std::unique_ptr<std::vector<size_t>> *indices)
		{
			auto new_indices = std::make_unique<std::vector<size_t>>(length);
			std::iota(new_indices->begin(), new_indices->end(), 0);

			auto val = std::stable_partition(new_indices->begin(), new_indices->end(), [begin, &flags](size_t index) { return (flags[index + begin] & 0x00FF) == 0; }) - new_indices->begin() + begin;

			*indices = std::move(new_indices);
			return val;
		}

		void HostParticleSnapshot::gather(const std::vector<size_t> &indices, size_t begin, size_t length)
		{
			sr::data::gather(r, indices, begin, length);
			sr::data::gather(v, indices, begin, length);
			sr::data::gather(id, indices, begin, length);
		}

		void HostParticleSnapshot::resize(size_t length)
		{
			r.resize(length);
			v.resize(length);
			id.resize(length);
			n_alive = std::min(length, n_alive);
			n = length;
		}

		void HostParticleSnapshot::filter(const std::vector<size_t> &filter, HostParticleSnapshot &out) const
		{
			out = HostParticleSnapshot(filter.size());
			out.n_alive = 0;

			size_t index = 0;

			for (size_t i : filter)
			{
				out.id[index] = id[i];
				out.r[index] = r[i];
				out.v[index] = v[i];

				if (i < n_alive)
				{
					out.n_alive++;
				}

				index++;
			}
		}

		void HostParticlePhaseSpace::gather(const std::vector<size_t> &indices, size_t begin, size_t length)
		{
			base.gather(indices, begin, length);
			sr::data::gather(deathtime(), indices, begin, length);
			sr::data::gather(deathflags(), indices, begin, length);

			if (deathtime_index().size() > 0)
			{
				sr::data::gather(deathtime_index(), indices, begin, length);
			}
		}

		void HostParticlePhaseSpace::filter(const std::vector<size_t> &filter, HostParticlePhaseSpace &out) const
		{
			base.filter(filter, out.base);
			out._deathflags = std::vector<uint16_t>(out.base.n);
			out._deathtime = std::vector<float>(out.base.n);

			size_t index = 0;

			for (size_t i : filter)
			{
				out._deathflags[index] = _deathflags[i];
				out._deathtime[index] = _deathtime[i];

				index++;
			}
		}

		std::unique_ptr<std::vector<size_t>> HostParticlePhaseSpace::stable_partition_unflagged(size_t begin, size_t length)
		{
			std::unique_ptr<std::vector<size_t>> indices;
			n_alive() = stable_partition_unflagged_indices(deathflags(), begin, length, &indices);

			this->gather(*indices, begin, length);

			return indices;
		}

		std::unique_ptr<std::vector<size_t>> HostParticlePhaseSpace::stable_partition_alive(size_t begin, size_t length)
		{
			std::unique_ptr<std::vector<size_t>> indices;
			n_alive() = stable_partition_alive_indices(deathflags(), begin, length, &indices);

			this->gather(*indices, begin, length);

			return indices;
		}

		std::unique_ptr<std::vector<size_t>> HostParticleSnapshot::sort_by_id(size_t begin, size_t length)
		{
			auto new_indices = std::make_unique<std::vector<size_t>>(length);

			std::iota(new_indices->begin(), new_indices->end(), 0);
			std::sort(new_indices->begin(), new_indices->end(), [begin, this](size_t a, size_t b) { return id[a + begin] < id[b + begin]; });

			this->gather(*new_indices, 0, length);

			return new_indices;
		}

		std::unique_ptr<std::vector<size_t>> HostParticlePhaseSpace::sort_by_id(size_t begin, size_t length)
		{
			auto indices = base.sort_by_id(begin, length);

			this->gather(*indices, 0, length);
			return indices;
		}

		Configuration::Configuration()
		{
			num_thread = 4;
			max_kep = 10;
			t_0 = 0;
			t_f = 365e4;
			dt = 122;
			tbsize = 1024;
			write_bary_track = false;
			keep_all_dumps = false;
			wh_ce_n1 = 8;
			wh_ce_n2 = 4;
			wh_ce_r1 = 1;
			wh_ce_r2 = 3.5;
			cull_radius = 0.5;

			resync_every = 1;
			print_every = 10;
			energy_every = 1;
			track_every = 0;
			split_track_file = 0;

			dump_every = 1000;
			max_particle = static_cast<uint32_t>(-1);
			big_g = 1;
			icsin = "";
			plin = "";
			hybridin = "";
			hybridout = "state.out";
			readsplit = 0;
			writesplit = 0;
			writebinary = 0;
			readbinary = 0;
			outfolder = "output/";
			readmomenta = false;
			writemomenta = false;
		}

		Configuration Configuration::output_config() const
		{
#pragma GCC warning "TODO"
			return *this;
		}

		void read_configuration(std::istream &in, Configuration *out)
		{
			size_t linenum = 0;
			std::string line;
			while (std::getline(in, line))
			{
				linenum++;
				if (line.length() == 0)
					continue;
				if (line[0] == '#')
					continue;

				size_t split = line.find(' ');
				if (split == std::string::npos)
				{
					std::ostringstream ss;
					ss << "Unrecognized line " << linenum << ": " << line;
					throw std::runtime_error(ss.str());
				}

				std::string first = line.substr(0, split);
				std::string second = line.substr(split + 1, line.length() - split - 1);

				try
				{
					if (first == "Initial-Time")
						out->t_0 = std::stod(second);
					else if (first == "Time-Step")
						out->dt = std::stod(second);
					else if (first == "Final-Time")
						out->t_f = std::stod(second);
					else if (first == "Time-Block-Size")
						out->tbsize = std::stou(second);
					else if (first == "Max-Kepler-Iterations")
						out->max_kep = std::stou(second);
					else if (first == "Big-G")
						out->big_g = std::stod(second);
					else if (first == "Cull-Radius")
						out->cull_radius = std::stod(second);
					else if (first == "CPU-Thread-Count")
						out->num_thread = std::stou(second);
					else if (first == "Limit-Particle-Count")
						out->max_particle = std::stou(second);
					else if (first == "Log-Interval")
						out->print_every = std::stou(second);
					else if (first == "Resync-Interval")
						out->resync_every = std::stou(second);
					else if (first == "Status-Interval")
						out->energy_every = std::stou(second);
					else if (first == "Track-Interval")
						out->track_every = std::stou(second);
					else if (first == "Write-Barycentric-Track")
						out->write_bary_track = std::stoi(second) != 0;
					else if (first == "Split-Track-File")
						out->split_track_file = std::stou(second);
					else if (first == "Dump-Interval")
						out->dump_every = std::stou(second);
					else if (first == "Keep-All-Dumps")
						out->keep_all_dumps = std::stoi(second) != 0;
					else if (first == "Write-Split-Output")
						out->writesplit = std::stoi(second) != 0;
					else if (first == "Read-Split-Input")
						out->readsplit = std::stoi(second) != 0;
					else if (first == "Write-Binary-Output")
						out->writebinary = std::stoi(second) != 0;
					else if (first == "Read-Binary-Input")
						out->readbinary = std::stoi(second) != 0;
					else if (first == "Input-File")
						out->hybridin = second;
					else if (first == "Output-Folder")
						out->outfolder = second;
					else if (first == "Particle-Input-File")
						out->icsin = second;
					else if (first == "Planet-Input-File")
						out->plin = second;
					else if (first == "Read-Input-Momenta")
						out->readmomenta = std::stoi(second) != 0;
					else if (first == "Write-Output-Momenta")
						out->writemomenta = std::stoi(second) != 0;
					else
						throw std::invalid_argument("bad");
				}
				catch (std::invalid_argument &)
				{
					std::ostringstream ss;
					ss << "Could not parse line " << linenum << ": " << line;
					throw std::runtime_error(ss.str());
				}
			}

			if (out->outfolder.empty())
			{
				std::cerr << "Warning: Output-Folder was not specified, defaulting to ./" << std::endl;
				out->outfolder = "./";
			}
			if (out->readsplit && out->readbinary)
			{
				std::cerr << "Warning: Read-Split-Input was selected but Read-Binary-Input was also specified. Ignoring Read-Binary-Input" << std::endl;
			}
			if (out->writesplit && out->writebinary)
			{
				std::cerr << "Warning: Write-Split-Input was selected but Write-Binary-Input was also specified. Ignoring Write-Binary-Input" << std::endl;
			}
			if (!out->writesplit && out->hybridin == "")
			{
				throw std::runtime_error("Error: Input-File was not specified");
			}
			if (out->readsplit && (out->plin == "" || out->icsin == ""))
			{
				throw std::runtime_error("Error: Read-Split-Input was selected but Particle-Input-File or Planet-Input-File were not specified");
			}
		}

		void write_configuration(std::ostream &outstream, const Configuration &out)
		{
			outstream << std::setprecision(17);
			outstream << "Initial-Time " << out.t_0 << std::endl;
			outstream << "Time-Step " << out.dt << std::endl;
			outstream << "Final-Time " << out.t_f << std::endl;
			outstream << "Time-Block-Size " << out.tbsize << std::endl;
			outstream << "Cull-Radius " << out.cull_radius << std::endl;
			outstream << "Max-Kepler-Iterations " << out.max_kep << std::endl;
			outstream << "CPU-Thread-Count " << out.num_thread << std::endl;
			outstream << "Limit-Particle-Count " << out.max_particle << std::endl;
			outstream << "Log-Interval " << out.print_every << std::endl;
			outstream << "Status-Interval " << out.energy_every << std::endl;
			outstream << "Track-Interval " << out.track_every << std::endl;
			outstream << "Resync-Interval " << out.resync_every << std::endl;
			outstream << "Write-Barycentric-Track " << out.write_bary_track << std::endl;
			outstream << "Split-Track-File " << out.split_track_file << std::endl;
			outstream << "Dump-Interval " << out.dump_every << std::endl;
			outstream << "Keep-All-Dumps " << out.keep_all_dumps << std::endl;
			outstream << "Write-Split-Output " << out.writesplit << std::endl;
			outstream << "Write-Binary-Output " << out.writebinary << std::endl;
			outstream << "Read-Split-Input " << out.readsplit << std::endl;
			outstream << "Read-Binary-Input " << out.readbinary << std::endl;
			outstream << "Input-File " << out.hybridin << std::endl;
			outstream << "Output-Folder " << out.outfolder << std::endl;
			outstream << "Particle-Input-File " << out.icsin << std::endl;
			outstream << "Planet-Input-File " << out.plin << std::endl;
			outstream << "Read-Input-Momenta " << out.readmomenta << std::endl;
			outstream << "Write-Output-Momenta " << out.writemomenta << std::endl;
		}

		bool load_planet_data(HostPlanetPhaseSpace &pl, const Configuration &config, std::istream &plin)
		{
			size_t npl;
			plin >> npl;

			pl = HostPlanetPhaseSpace(npl, config.tbsize);

			for (size_t i = 0; i < npl; i++)
			{
				plin >> pl.m()[i];
				plin >> pl.r()[i].x >> pl.r()[i].y >> pl.r()[i].z;
				plin >> pl.v()[i].x >> pl.v()[i].y >> pl.v()[i].z;

				pl.id()[i] = static_cast<uint32_t>(i);
			}

			return false;
		}

		bool load_data_nohybrid(HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, const Configuration &config, std::istream &plin, std::istream &icsin)
		{
			load_planet_data(pl, config, plin);

			size_t npart;
			icsin >> npart;
			npart = std::min(npart, static_cast<size_t>(config.max_particle));

			pa = HostParticlePhaseSpace(npart);

			for (size_t i = 0; i < npart; i++)
			{
				icsin >> pa.r()[i].x >> pa.r()[i].y >> pa.r()[i].z;
				icsin >> pa.v()[i].x >> pa.v()[i].y >> pa.v()[i].z;

				std::string s;
				icsin >> s;
				if (!isdigit(s[0]))
				{
					icsin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					pa.deathtime()[i] = 0;
					pa.id()[i] = static_cast<uint32_t>(i);
					pa.deathflags()[i] = 0;
				}
				else
				{
					pa.deathtime()[i] = std::stof(s);
					icsin >> pa.deathflags()[i] >> pa.id()[i];
				}
			}

			return false;
		}

		bool load_data_hybrid(HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, const Configuration &config, std::istream &in)
		{
			std::string s;

			std::getline(in, s);
			std::istringstream ss(s);
			size_t npl;
			ss >> npl;

			pl = HostPlanetPhaseSpace(npl, config.tbsize);

			for (size_t i = 0; i < npl; i++)
			{
				std::getline(in, s);
				ss = std::istringstream(s);
				ss >> pl.m()[i];

				std::getline(in, s);
				ss = std::istringstream(s);
				ss >> pl.r()[i].x >> pl.r()[i].y >> pl.r()[i].z;

				std::getline(in, s);
				ss = std::istringstream(s);
				ss >> pl.v()[i].x >> pl.v()[i].y >> pl.v()[i].z;

				std::getline(in, s);
				ss = std::istringstream(s);

				ss >> pl.id()[i];
			}

			std::getline(in, s);
			ss = std::istringstream(s);
			size_t npart;
			ss >> npart;
			npart = std::min(npart, static_cast<size_t>(config.max_particle));
			pa = HostParticlePhaseSpace(npart);

			for (size_t i = 0; i < npart; i++)
			{
				std::getline(in, s);
				ss = std::istringstream(s);
				ss >> pa.r()[i].x >> pa.r()[i].y >> pa.r()[i].z;

				std::getline(in, s);
				ss = std::istringstream(s);
				ss >> pa.v()[i].x >> pa.v()[i].y >> pa.v()[i].z;

				std::getline(in, s);
				ss = std::istringstream(s);
				ss >> pa.id()[i] >> pa.deathflags()[i] >> pa.deathtime()[i];
			}

			return false;
		}

		bool load_data_hybrid_binary(HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, const Configuration &config, std::istream &in)
		{
			uint64_t templl;
			read_binary<uint64_t>(in, templl);
			size_t npl = static_cast<size_t>(templl);
			pl = HostPlanetPhaseSpace(npl, config.tbsize);

			for (size_t i = 0; i < pl.n(); i++)
			{
				read_binary<uint32_t>(in, pl.id()[i]);
				read_binary<double>(in, pl.m()[i]);
				read_binary<double>(in, pl.r()[i].x);
				read_binary<double>(in, pl.r()[i].y);
				read_binary<double>(in, pl.r()[i].z);
				read_binary<double>(in, pl.v()[i].x);
				read_binary<double>(in, pl.v()[i].y);
				read_binary<double>(in, pl.v()[i].z);
			}

			read_binary<uint64_t>(in, templl);
			size_t npart = std::min(static_cast<size_t>(templl), static_cast<size_t>(config.max_particle));
			pa = HostParticlePhaseSpace(npart);

			for (size_t i = 0; i < pa.n(); i++)
			{
				read_binary<uint32_t>(in, pa.id()[i]);
				read_binary<double>(in, pa.r()[i].x);
				read_binary<double>(in, pa.r()[i].y);
				read_binary<double>(in, pa.r()[i].z);
				read_binary<double>(in, pa.v()[i].x);
				read_binary<double>(in, pa.v()[i].y);
				read_binary<double>(in, pa.v()[i].z);
				read_binary<uint16_t>(in, pa.deathflags()[i]);
				read_binary<float>(in, pa.deathtime()[i]);
			}

			return !in;
		}

		bool load_data(HostPlanetPhaseSpace &pl, HostParticlePhaseSpace &pa, const Configuration &config)
		{
			bool ret;
			if (config.readsplit)
			{
				if (!sr::util::does_file_exist(config.plin))
				{
					std::ostringstream ss;
					ss << "Planet input file " << config.plin << " does not exist";
					throw std::runtime_error(ss.str());
				}
				if (!sr::util::does_file_exist(config.icsin))
				{
					std::ostringstream ss;
					ss << "Particle input file " << config.icsin << " does not exist";
					throw std::runtime_error(ss.str());
				}

				std::ifstream plinfile(config.plin), icsinfile(config.icsin);
				ret = load_data_nohybrid(pl, pa, config, plinfile, icsinfile);
			}
			else
			{
				if (!sr::util::does_file_exist(config.hybridin))
				{
					std::ostringstream ss;
					ss << "Input file " << config.hybridin << " does not exist";
					throw std::runtime_error(ss.str());
				}

				if (config.readbinary)
				{
					std::ifstream in(config.hybridin, std::ios_base::binary);
					ret = load_data_hybrid_binary(pl, pa, config, in);
				}
				else
				{
					std::ifstream in(config.hybridin);
					ret = load_data_hybrid(pl, pa, config, in);
				}
			}

			if (!ret)
			{
				for (size_t i = 0; i < pl.n(); i++)
				{
					if (config.readmomenta)
					{
						pl.v()[i] /= pl.m()[i];
					}

					pl.m()[i] *= config.big_g;
				}

				pa.stable_partition_alive(0, pa.n());
			}

			return ret;
		}

		void save_data_hybrid_binary(const HostPlanetSnapshot &pl, const HostParticlePhaseSpace &pa, const Configuration &config, std::ostream &out)
		{
			(void)config;
			write_binary(out, static_cast<uint64_t>(pl.n_alive));
			for (size_t i = 0; i < pl.n_alive; i++)
			{
				double m = pl.m[i];
				write_binary(out, pl.id[i]);
				write_binary(out, m);

				if (!config.writemomenta)
					m = 1;
				write_binary(out, pl.r[i].x);
				write_binary(out, pl.r[i].y);
				write_binary(out, pl.r[i].z);
				write_binary(out, pl.v[i].x * m);
				write_binary(out, pl.v[i].y * m);
				write_binary(out, pl.v[i].z * m);
			}

			write_binary(out, static_cast<uint64_t>(pa.n()));
			for (size_t i = 0; i < pa.n(); i++)
			{
				write_binary(out, pa.id()[i]);
				write_binary(out, pa.r()[i].x);
				write_binary(out, pa.r()[i].y);
				write_binary(out, pa.r()[i].z);
				write_binary(out, pa.v()[i].x);
				write_binary(out, pa.v()[i].y);
				write_binary(out, pa.v()[i].z);
				write_binary(out, pa.deathflags()[i]);
				write_binary(out, pa.deathtime()[i]);
			}
		}

		void save_data_hybrid(const HostPlanetSnapshot &pl, const HostParticlePhaseSpace &pa, const Configuration &config, std::ostream &out)
		{
			(void)config;
			out << pl.n_alive << std::endl;
			out << std::setprecision(17);
			for (size_t i = 0; i < pl.n_alive; i++)
			{
				double m = pl.m[i];
				out << m << std::endl;
				if (!config.writemomenta)
					m = 1;
				out << pl.r[i].x << " " << pl.r[i].y << " " << pl.r[i].z << std::endl;
				out << pl.v[i].x * m << " " << pl.v[i].y * m << " " << pl.v[i].z * m << std::endl;
				out << pl.id[i] << std::endl;
			}

			out << pa.n() << std::endl;
			out << std::setprecision(17);
			for (size_t i = 0; i < pa.n(); i++)
			{
				out << pa.r()[i].x << " " << pa.r()[i].y << " " << pa.r()[i].z << std::endl;
				out << pa.v()[i].x << " " << pa.v()[i].y << " " << pa.v()[i].z << std::endl;
				out << pa.id()[i] << " " << pa.deathflags()[i] << " " << pa.deathtime()[i] << std::endl;
			}
		}

		void save_data_nohybrid(const HostPlanetSnapshot &pl, const HostParticlePhaseSpace &pa, const Configuration &config, std::ostream &plout, std::ostream &icsout)
		{
			(void)config;
			plout << pl.n_alive << std::endl;
			plout << std::setprecision(17);
			for (size_t i = 0; i < pl.n_alive; i++)
			{
				double m = pl.m[i];
				plout << m << std::endl;

				if (!config.writemomenta)
					m = 1;
				plout << pl.r[i].x << " " << pl.r[i].y << " " << pl.r[i].z << std::endl;
				plout << pl.v[i].x * m << " " << pl.v[i].y * m << " " << pl.v[i].z * m << std::endl;
			}

			icsout << pa.n() << std::endl;
			icsout << std::setprecision(17);
			for (size_t i = 0; i < pa.n(); i++)
			{
				icsout << pa.r()[i].x << " " << pa.r()[i].y << " " << pa.r()[i].z << std::endl;
				icsout << pa.v()[i].x << " " << pa.v()[i].y << " " << pa.v()[i].z << std::endl;
				icsout << pa.deathtime()[i] << " " << pa.deathflags()[i] << " " << pa.id()[i] << std::endl;
			}
		}

		void save_data_swift(const HostPlanetSnapshot &pl, const HostParticlePhaseSpace &pa, std::ostream &plout, std::ostream &icsout)
		{
			plout << pl.n_alive << std::endl;
			plout << std::setprecision(17);
			for (size_t i = 0; i < pl.n_alive; i++)
			{
				double m = pl.m[i];
				plout << m << std::endl;

				plout << pl.r[i].x << " " << pl.r[i].y << " " << pl.r[i].z << std::endl;
				plout << pl.v[i].x << " " << pl.v[i].y << " " << pl.v[i].z << std::endl;
			}

			icsout << pa.n() << std::endl;
			icsout << std::setprecision(17);
			for (size_t i = 0; i < pa.n(); i++)
			{
				icsout << pa.r()[i].x << " " << pa.r()[i].y << " " << pa.r()[i].z << std::endl;
				icsout << pa.v()[i].x << " " << pa.v()[i].y << " " << pa.v()[i].z << std::endl;
				icsout << "0 0 0 0 0 0 0 0 0 0 0 0 0" << std::endl;
				icsout << "0.0 0.0 0.0 0.0 0.0" << std::endl
					   << "0.0 0.0 0.0 0.0 0.0" << std::endl
					   << "0.0 0.0 0.0" << std::endl;
			}
		}

		void save_data(const HostPlanetSnapshot &pl, const HostParticlePhaseSpace &pa, const Configuration &config, const std::string &outfile)
		{
			if (!config.writesplit)
			{
				std::ostringstream ss;

				if (config.writebinary)
				{
					std::ofstream out(outfile, std::ios_base::binary);
					save_data_hybrid_binary(pl, pa, config, out);
				}
				else
				{
					std::ofstream out(outfile);
					save_data_hybrid(pl, pa, config, out);
				}
			}
			else
			{
				std::ofstream ploutfile(sr::util::joinpath(config.outfolder, "pl.out")), icsoutfile(sr::util::joinpath(config.outfolder, "ics.out"));
				save_data_nohybrid(pl, pa, config, ploutfile, icsoutfile);
			}
		}

		void save_binary_track(std::ostream &trackout, const HostPlanetSnapshot &pl, const HostParticleSnapshot &pa, double time, bool to_elements, bool barycentric_elements)
		{
			sr::data::write_binary(trackout, static_cast<double>(time));

			if (pl.n_alive > 0)
			{
				sr::data::write_binary(trackout, static_cast<uint64_t>(pl.n_alive));
			}
			else
			{
				sr::data::write_binary(trackout, static_cast<uint64_t>(0));
			}

			for (uint32_t i = 0; i < pl.n_alive; i++)
			{
				if (to_elements)
				{
					double center_mass = pl.m[0];
					f64_3 center_r = pl.r[0] * pl.m[0];
					f64_3 center_v = pl.v[0] * pl.m[0];
					if (barycentric_elements)
					{
						for (uint32_t j = 1; j < pl.n_alive; j++)
						{
							if (j != i)
							{
								center_mass += pl.m[j];
								center_r += pl.r[j] * pl.m[j];
								center_v += pl.v[j] * pl.m[j];
							}
						}
					}
					center_r /= center_mass;
					center_v /= center_mass;

					double a, e, in, capom, om, f;
					sr::convert::to_elements(pl.m[i] + center_mass, pl.r[i] - center_r, pl.v[i] - center_v,
											 nullptr, &a, &e, &in, &capom, &om, &f);

					sr::data::write_binary(trackout, static_cast<uint32_t>(pl.id[i]));
					sr::data::write_binary(trackout, static_cast<float>(a));
					sr::data::write_binary(trackout, static_cast<float>(e));
					sr::data::write_binary(trackout, static_cast<float>(in));
					sr::data::write_binary(trackout, static_cast<float>(capom));
					sr::data::write_binary(trackout, static_cast<float>(om));
					sr::data::write_binary(trackout, static_cast<float>(f));
				}
				else
				{
					sr::data::write_binary(trackout, static_cast<uint32_t>(pl.id[i]));
					sr::data::write_binary(trackout, static_cast<float>(pl.r[i].x));
					sr::data::write_binary(trackout, static_cast<float>(pl.r[i].y));
					sr::data::write_binary(trackout, static_cast<float>(pl.r[i].z));
					sr::data::write_binary(trackout, static_cast<float>(pl.v[i].x));
					sr::data::write_binary(trackout, static_cast<float>(pl.v[i].y));
					sr::data::write_binary(trackout, static_cast<float>(pl.v[i].z));
				}
			}

			sr::data::write_binary(trackout, static_cast<uint64_t>(pa.n_alive));
			for (uint32_t i = 0; i < pa.n_alive; i++)
			{
				if (to_elements)
				{
					double center_mass = pl.m[0];
					f64_3 center_r = pl.r[0];
					f64_3 center_v = pl.v[0];
					if (barycentric_elements)
					{
						sr::convert::find_barycenter(pl.r, pl.v, pl.m, pl.n_alive, center_r, center_v, center_mass);
					}

					double a, e, in, capom, om, f;
					sr::convert::to_elements(center_mass, pa.r[i] - center_r, pa.v[i] - center_v,
											 nullptr, &a, &e, &in, &capom, &om, &f);

					sr::data::write_binary(trackout, static_cast<uint32_t>(pa.id[i]));
					sr::data::write_binary(trackout, static_cast<float>(a));
					sr::data::write_binary(trackout, static_cast<float>(e));
					sr::data::write_binary(trackout, static_cast<float>(in));
					sr::data::write_binary(trackout, static_cast<float>(capom));
					sr::data::write_binary(trackout, static_cast<float>(om));
					sr::data::write_binary(trackout, static_cast<float>(f));
				}
				else
				{
					sr::data::write_binary(trackout, static_cast<uint32_t>(pa.id[i]));
					sr::data::write_binary(trackout, static_cast<float>(pa.r[i].x));
					sr::data::write_binary(trackout, static_cast<float>(pa.r[i].y));
					sr::data::write_binary(trackout, static_cast<float>(pa.r[i].z));
					sr::data::write_binary(trackout, static_cast<float>(pa.v[i].x));
					sr::data::write_binary(trackout, static_cast<float>(pa.v[i].y));
					sr::data::write_binary(trackout, static_cast<float>(pa.v[i].z));
				}
			}

			trackout.flush();
		}

		void TrackReader::check_state(const State &expected)
		{
			if (state != expected)
			{
				throw std::runtime_error("invalid reader state");
			}
		}

		void TrackReader::read_time()
		{
			check_state(State::Start);

			sr::data::read_binary<double>(input, time);
			state = State::PlanetsBegin;
		}

		void TrackReader::begin_planets()
		{
			check_state(State::PlanetsBegin);

			planets = HostPlanetSnapshot();

			uint64_t templl;
			sr::data::read_binary<uint64_t>(input, templl);

			n_planets = static_cast<size_t>(templl);
			state = State::PlanetsEnd;
		}

		void TrackReader::read_planets(const std::vector<uint32_t> *filter)
		{
			check_state(State::PlanetsEnd);
			std::streampos pos = input.tellg();

			// planets = HostPlanetSnapshot(static_cast<size_t>(n_planets));
			planets = HostPlanetSnapshot(static_cast<size_t>(filter ? filter->size() : n_planets));

			int index = 0;
			for (uint32_t i = 0; i < n_planets; i++)
			{
				uint32_t id;
				sr::data::read_binary<uint32_t>(input, id);

				if (!filter || std::find(filter->begin(), filter->end(), id) != filter->end())
				{
					planets.id[index] = id;

					float tempfloat;
					sr::data::read_binary<float>(input, tempfloat);
					planets.r[index].x = tempfloat;
					sr::data::read_binary<float>(input, tempfloat);
					planets.r[index].y = tempfloat;
					sr::data::read_binary<float>(input, tempfloat);
					planets.r[index].z = tempfloat;
					sr::data::read_binary<float>(input, tempfloat);
					planets.v[index].x = tempfloat;
					sr::data::read_binary<float>(input, tempfloat);
					planets.v[index].y = tempfloat;
					sr::data::read_binary<float>(input, tempfloat);
					planets.v[index].z = tempfloat;
					index++;
				}
				else
				{
					input.seekg(4 * 6, std::ios_base::cur);
				}
			}

			input.seekg(pos);
		}

		void TrackReader::end_planets()
		{
			check_state(State::PlanetsEnd);

			input.seekg(TRACK_PLANET_STRIDE * n_planets, std::ios_base::cur);
			state = State::ParticlesBegin;
		}

		void TrackReader::begin_particles()
		{
			check_state(State::ParticlesBegin);

			particles = HostParticleSnapshot();

			uint64_t templl;
			sr::data::read_binary<uint64_t>(input, templl);

			n_particles = static_cast<size_t>(templl);
			state = State::ParticlesEnd;
		}

		void TrackReader::read_particles()
		{
			check_state(State::ParticlesEnd);

			particles = HostParticleSnapshot(n_particles);
			std::streampos position = input.tellg();

			for (uint32_t i = 0; i < n_particles; i++)
			{
				sr::data::read_binary<uint32_t>(input, particles.id[i]);

				float tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				particles.r[i].x = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				particles.r[i].y = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				particles.r[i].z = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				particles.v[i].x = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				particles.v[i].y = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				particles.v[i].z = tempfloat;
			}

			input.seekg(position);
		}

		bool TrackReader::read_particle(uint32_t id)
		{
			check_state(State::ParticlesEnd);

			std::streampos position = input.tellg();

			size_t result = bsearch_track(input, n_particles, id, TRACK_PARTICLE_STRIDE);

			if (result)
			{
				size_t newsize = particles.n + 1;
				particles.n = newsize;
				particles.n_alive = newsize;

				particles.r.resize(newsize);
				particles.v.resize(newsize);
				particles.id.resize(newsize);

				particles.r[newsize - 1].x = sr::data::read_binary<float>(input);
				particles.r[newsize - 1].y = sr::data::read_binary<float>(input);
				particles.r[newsize - 1].z = sr::data::read_binary<float>(input);
				particles.v[newsize - 1].x = sr::data::read_binary<float>(input);
				particles.v[newsize - 1].y = sr::data::read_binary<float>(input);
				particles.v[newsize - 1].z = sr::data::read_binary<float>(input);

				particles.id[newsize - 1] = id;
			}

			input.seekg(position);

			return result > 0;
		}

		void TrackReader::end_particles()
		{
			check_state(State::ParticlesEnd);
			input.seekg(TRACK_PARTICLE_STRIDE * n_particles, std::ios_base::cur);

			state = State::Finish;
		}

		size_t TrackReader::bsearch_track(std::istream &f, size_t npa, uint32_t partnum, size_t stride)
		{
			int64_t base = f.tellg();
			int64_t left = 0;
			int64_t right = npa - 1;

			size_t its = 0;
			while (left <= right)
			{
				its++;
				int64_t mid = left + (right - left) / 2;
				f.seekg(base + mid * stride);

				uint32_t id;
				sr::data::read_binary<uint32_t>(f, id);
				if (id == partnum)
				{
					return its;
				}
				else if (id > partnum)
				{
					right = mid - 1;
				}
				else
				{
					left = mid + 1;
				}
			}

			return 0;
		}

		void process_track(std::istream &input,
						   const TrackReaderOptions &options,
						   const std::function<void(HostPlanetSnapshot &, HostParticleSnapshot &, double)> &callback)
		{
			while (true)
			{
				TrackReader reader(input);
				bool skip = false;

				reader.read_time();
				if (!input)
					break;

				if (reader.time > options.max_time)
				{
					skip = true;
				}

				reader.begin_planets();
				if (!input)
					break;

				if (!skip)
				{
					reader.read_planets(options.take_all_planets ? nullptr : &options.planet_filter);
				}
				reader.end_planets();

				reader.begin_particles();
				if (!input)
					break;

				if (!skip)
				{
					if (options.take_all_particles)
					{
						reader.read_particles();
					}
					else
					{
						for (uint32_t id : options.particle_filter)
						{
							reader.read_particle(id);
						}
					}
				}
				reader.end_particles();

				if (!input)
					break;

				if (!skip)
				{
					callback(reader.planets, reader.particles, reader.time);
				}
			}
		}

		void read_tracks(const std::string &path,
						 const TrackReaderOptions &options,
						 const std::function<void(HostPlanetSnapshot &, HostParticleSnapshot &, double)> &callback)
		{
			sr::util::PathType pathtype = sr::util::get_path_type(path);

			if (pathtype == sr::util::PathType::Directory)
			{
				for (size_t i = 0; true; i++)
				{
					std::ostringstream ss;
					ss << path;
					if (path[path.size() - 1] != '/')
						ss << '/';
					ss << "track." << i << ".out";

					if (!options.silent)
					{
						std::cout << "Reading " << ss.str() << std::endl;
					}

					if (!sr::util::does_file_exist(ss.str()))
					{
						if (i == 0)
						{
							throw std::runtime_error("Directory not contain track files");
						}
						else
						{
							if (!options.silent)
							{
								std::cout << i << " files read" << std::endl;
							}
							return;
						}
					}

					std::ifstream input(ss.str(), std::ios_base::binary);

					process_track(input, options, callback);
				}
			}
			else if (pathtype == sr::util::PathType::File)
			{
				std::ifstream input(path, std::ios_base::binary);

				process_track(input, options, callback);
			}
			else
			{
				throw std::runtime_error("Path does not exist");
			}
		}
	} // namespace data
} // namespace sr
