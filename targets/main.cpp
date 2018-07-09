/*************************************************************

.---. .            .         .   .-,--.                        
\___  |- ,-. ,-. ,-| . . ,-. |-   `|__/ ,-. .  , ,-. ,-. . ,-. 
    \ |  ,-| |   | | | | `-. |    /  \  |-' | /  |-' |   | |-' 
`---' `' `-^ '   `-^ `-^ `-' `'   `-  ` `-' `'   `-' '   ' `-' 

*************************************************************/

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <cmath>
#include <iomanip>

#include <execinfo.h>
#include <sys/stat.h>
#include <csignal>

#ifdef NO_CUDA
	#include "../cpu_executor.h"
#else
	#include "../executor_facade.h"
#endif

#include "../data.h"
#include "../wh.h"
#include "../convert.h"
#include "../util.h"


volatile sig_atomic_t end_loop = 0;

void term(int signum)
{
	(void) signum;
	end_loop = 1;
}

int main(int argv, char** argc)
{
	std::ios_base::sync_with_stdio(false);

	if (sizeof(double) != 8)
	{
		std::cerr << "sizeof double is not 8 - this code will likely not work!" << std::endl;
	}

	if (sizeof(float) != 4)
	{
		std::cerr << "sizeof float is not 4 - this code will likely not work!" << std::endl;
	}

	std::string configin = "config.in";
	if (argv >= 2) configin = std::string(argc[1]);
	
	std::cout << "Reading from configuration file " << configin << std::endl;
	
	std::ifstream configfile(configin);

	sr::data::Configuration config_mut;
	if (read_configuration(configfile, &config_mut)) return -1;

	const sr::data::Configuration& config = config_mut;

	{
		std::ofstream configstream(sr::util::joinpath(config.outfolder, "config.in"));
		write_configuration(configstream, config);
	}


#ifdef NO_CUDA
	if (config_mut.use_gpu)
	{
		std::cout << "Enable-GPU cannot be enabled when compiling without CUDA!" << std::endl;
		return -1;
	}
#endif


	sr::data::Configuration out_config = config.output_config();

	mkdir(config.outfolder.c_str(), ACCESSPERMS);

	if (!sr::util::is_dir_empty(config.outfolder))
	{
		std::cout << "Output folder is not empty! Do you want to continue?" << std::endl;
		std::cout << "Type \"Yes\" exactly as shown to continue: ";
	
		std::string s;
		std::getline(std::cin, s);

		if (s != "Yes") return -1;
	}

	mkdir(sr::util::joinpath(config.outfolder, "dump").c_str(), ACCESSPERMS);
	mkdir(sr::util::joinpath(config.outfolder, "tracks").c_str(), ACCESSPERMS);

	std::ofstream coutlog(sr::util::joinpath(config.outfolder, "stdout"));
	sr::util::teestream tout(std::cout, coutlog);

	tout << "Host uses little-endian floats? " << (sr::data::is_float_little_endian() ? "yes" : "no") << std::endl;
	tout << "Host uses little-endian doubles? " << (sr::data::is_double_little_endian() ? "yes" : "no") << std::endl;
	tout << "Host uses little-endian ints? " << (sr::data::is_int_little_endian() ? "yes" : "no") << std::endl;

	sr::data::HostData hd;

#ifdef NO_CUDA
	sr::exec::CPUExecutor ex(hd, config, tout);
#else
	sr::exec::ExecutorFacade ex(hd, config, tout);
#endif

	ex.t = config.t_0;

	if (load_data(hd, config)) return -1;

	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);

	std::ofstream encounterlog(sr::util::joinpath(config.outfolder, "encounter.out"));
	ex.encounter_output = &encounterlog;

	std::ofstream timelog(sr::util::joinpath(config.outfolder, "time.out"));
	timelog << "start " << std::put_time(&tm, "%c %Z") << std::endl;

	ex.init();

	uint32_t counter = 0;
	uint32_t dump_num = 0;

	bool crashed = false;
	std::ofstream trackout;

	signal(SIGTERM, term);
	signal(SIGINT, term);

	uint32_t track_num = 1;

	try
	{
		trackout = std::ofstream(sr::util::joinpath(config.outfolder, "tracks/track.0.out"), std::ios_base::binary);

		while (ex.t < config.t_f)
		{
#ifdef NO_CUDA
			double cputimeout;
		       	ex.loop(&cputimeout);

			double timediff = cputimeout;
#else
			double cputimeout, gputimeout;
		       	ex.loop(&cputimeout, &gputimeout);

			double timediff = gputimeout - cputimeout;
#endif

			counter++;

			ex.add_job([&timelog, &tout, &ex, &config, counter, timediff]()
				{
					bool output_energy = config.energy_every != 0 && (counter % config.energy_every == 0);
					bool log_out = config.print_every != 0 && (counter % config.print_every == 0);


					if (!log_out && !output_energy) return;

					double e_;
					f64_3 l_;
					sr::wh::calculate_planet_metrics(ex.hd.planets, &e_, &l_);
					double elapsed = ex.time();
					double total = elapsed * (config.t_f - config.t_0) / (ex.t - config.t_0);

					if (output_energy)
					{
						timelog << std::setprecision(13) << "time " << elapsed << " " << ex.hd.particles.n_alive << " " << ex.hd.particles.n_encounter << std::endl;
						timelog << "ep " << e_ << std::endl;
						timelog << "lp " << l_ << std::endl;
					}

					if (log_out)
					{
						tout << std::setprecision(4);
						tout << "t=" << ex.t << " (" << elapsed / total * 100 << "% " << elapsed << "m elapsed, "
							<< total << "m total " << total - elapsed << "m remain)" << std::endl;
						tout << "Error = " << (e_ - ex.e_0) / ex.e_0 * 100 << ", " <<
							ex.hd.particles.n_alive << " particles remaining, " << ex.hd.particles.n_encounter << " in encounter" << std::endl;

						tout << "GPU took " << std::setprecision(4) << timediff << " ms longer than CPU" << std::endl;
					}
				});
			
			bool dump = config.dump_every != 0 && counter % config.dump_every == 0;
			bool track = config.track_every != 0 && counter % config.track_every == 0;

			if (dump || track)
			{
#ifndef NO_CUDA
				ex.download_data();
#endif

				if (dump)
				{
					out_config.t_f = config.t_f - config.t_0 + ex.t;
					out_config.t_0 = ex.t;
					ex.add_job([&tout, &ex, &out_config, &config, &dump_num]()
						{
							tout << "Dumping to disk. t = " << ex.t << std::endl;
							std::ostringstream ss;
							ss << "dumps/config." << dump_num << ".out";

							std::ofstream configout(sr::util::joinpath(config.outfolder, ss.str()));
							write_configuration(configout, out_config);

							ss = std::ostringstream();
							ss << "dumps/state." << dump_num << ".out";
							save_data(ex.hd, config, sr::util::joinpath(config.outfolder, ss.str()), true);

							dump_num++;
						});
				}

				if (track)
				{
					if (config.split_track_file > 0 && trackout.tellp() > config.split_track_file)
					{
						std::ostringstream ss;
						ss << "tracks/track." << track_num++ << ".out";
						trackout = std::ofstream(sr::util::joinpath(config.outfolder, ss.str()), std::ios_base::binary);
					}

					ex.add_job([&trackout, &ex, &config]()
						{
							sr::data::write_binary(trackout, ex.t);
							sr::data::write_binary(trackout, ex.hd.planets.n_alive - 1);

							for (uint32_t i = 1; i < ex.hd.planets.n_alive; i++)
							{
								double a, e, in, capom, om, f;
								sr::convert::to_elements(ex.hd.planets.m[i] + ex.hd.planets.m[0], ex.hd.planets.r[i], ex.hd.planets.v[i],
									nullptr, &a, &e, &in, &capom, &om, &f);

								sr::data::write_binary(trackout, static_cast<uint32_t>(ex.hd.planets.id[i]));
								sr::data::write_binary(trackout, static_cast<float>(a));
								sr::data::write_binary(trackout, static_cast<float>(e));
								sr::data::write_binary(trackout, static_cast<float>(in));
								sr::data::write_binary(trackout, static_cast<float>(capom));
								sr::data::write_binary(trackout, static_cast<float>(om));
								sr::data::write_binary(trackout, static_cast<float>(f));
							}

							sr::data::write_binary(trackout, ex.hd.particles.n_alive);
							for (uint32_t i = 0; i < ex.hd.particles.n_alive; i++)
							{
								double a, e, in, capom, om, f;
								sr::convert::to_elements(ex.hd.planets.m[0], ex.hd.particles.r[i], ex.hd.particles.v[i],
									nullptr, &a, &e, &in, &capom, &om, &f);

								sr::data::write_binary(trackout, static_cast<uint32_t>(ex.hd.particles.id[i]));
								sr::data::write_binary(trackout, static_cast<float>(a));
								sr::data::write_binary(trackout, static_cast<float>(e));
								sr::data::write_binary(trackout, static_cast<float>(in));
								sr::data::write_binary(trackout, static_cast<float>(capom));
								sr::data::write_binary(trackout, static_cast<float>(om));
								sr::data::write_binary(trackout, static_cast<float>(f));
							}

							trackout.flush();
						});
				}
			}


			if (end_loop)
			{
				tout << "Caught signal." << std::endl;
				throw std::exception();
			}
		}
	}
	catch (const std::exception& e)
	{
		void* array[50];
		size_t size = backtrace(array, 50);
		backtrace_symbols_fd(array, static_cast<int>(size), 1);

		tout << "Exception caught: " << std::endl;
		tout << e.what() << std::endl;
		tout << "Recovering data." << std::endl;
		crashed = true;
	}

	ex.finish();
#ifndef NO_CUDA
	ex.download_data(true);
#endif

	tout << "Saving to disk." << std::endl;
	save_data(hd, config, sr::util::joinpath(config.outfolder, "state.out"));
	out_config.t_f = config.t_f - config.t_0 + ex.t;
	out_config.t_0 = ex.t;
	std::ofstream configout(sr::util::joinpath(config.outfolder, "config.out"));
	write_configuration(configout, out_config);

	t = std::time(nullptr);
	tm = *std::localtime(&t);
	timelog << "end " << std::put_time(&tm, "%c %Z") << std::endl;

	return crashed ? -1 : 0;
}
