#include "wh.h"
#include <ctime>
#include <chrono>
#include <functional>
#include <ostream>

namespace sr
{
namespace exec
{
	struct CPUExecutor
	{
		sr::data::HostData& hd;
		sr::wh::WHIntegrator integrator;

		float64_t t;
		float64_t e_0;

		std::ostream& output;

		const sr::data::Configuration& config;

		std::chrono::time_point<std::chrono::high_resolution_clock> starttime;

		std::vector<std::function<void()>> work;

		CPUExecutor(const CPUExecutor&) = delete;
		CPUExecutor(sr::data::HostData& hd, const sr::data::Configuration& config, std::ostream& out);

		void init();

		double time() const;
		void loop(double* cputime);
		void add_job(const std::function<void()>& job);

		void resync();
		void finish();
		void swap_logs();
		void step_planets();
	};
}
}
