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
#include <csignal>

#ifdef NO_CUDA
	#include "../src/cpu_executor.h"
#else
	#include "../src/executor_facade.h"
#endif

#include "../src/data.h"
#include "../src/wh.h"
#include "../src/convert.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

static const char USAGE[] = R"(sr(_cpu)
Usage:
    sr(_cpu) [options] [<config>]

Options:
    -h, --help         Show this screen.
)";

volatile sig_atomic_t end_loop = 0;

void term(int signum)
{
	(void) signum;
	end_loop = 1;
}

int main(int argc, char** argv)
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

	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "sr");

	std::string configin = "config.in";
	if (args["<config>"]) configin = args["<config>"].asString();
	
	std::cout << "Reading from configuration file " << configin << std::endl;
	
	std::ifstream configfile(configin);

	sr::data::Configuration config_mut;

	try
	{
		read_configuration(configfile, &config_mut);
	}
	catch (std::exception& e)
	{
		std::cerr << "Could not read configuration." << std::endl;
		std::cerr << e.what() << std::endl;
		return -1;
	}

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

	sr::util::make_dir(config.outfolder);

	if (!sr::util::is_dir_empty(config.outfolder))
	{
		std::cout << "Output folder is not empty! Do you want to continue?" << std::endl;
		std::cout << "Type \"Yes\" exactly as shown to continue: ";
	
		std::string s;
		std::getline(std::cin, s);

		if (s != "Yes") return -1;
	}

	sr::util::make_dir(sr::util::joinpath(config.outfolder, "dump"));
	sr::util::make_dir(sr::util::joinpath(config.outfolder, "tracks"));

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
							sr::data::save_binary_track(trackout, ex.hd.planets_snapshot, ex.hd.particles.make_snapshot(), ex.t, true);
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
