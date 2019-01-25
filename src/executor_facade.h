#pragma once
#include "data.h"
#include "wh.h"
#include <functional>
#include <memory>

namespace sr
{
	namespace data
	{
		struct DeviceData;
	}

	namespace exec
	{
		struct Executor;

		struct ExecutorFacade
		{
			sr::data::HostData& hd;
			std::unique_ptr<sr::data::DeviceData> dd;
			std::unique_ptr<Executor> impl;

			float64_t& t;
			float64_t& e_0;


			ExecutorFacade(sr::data::HostData& hd, const sr::data::Configuration& config, std::ostream& out);
			~ExecutorFacade();

			void init();
			void download_data(bool ignore_errors = false);
			double time() const;
			void loop(double* cputimeout, double* gputimeout);
			void add_job(const std::function<void()>& job);
			void finish();
		};
	}
}
