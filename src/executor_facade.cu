#include "executor_facade.h"
#include "util.h"
#include "data.cuh"
#include "executor.cuh"

namespace sr
{
namespace exec
{
	using namespace sr::data;

	ExecutorFacade::ExecutorFacade(HostData& _hd, const Configuration& config, std::ostream& out) :
		hd(_hd),
		dd(std::make_unique<DeviceData>()),
		impl(std::make_unique<Executor>(_hd, *dd.get(), config, out)),
		t(impl->t),
		e_0(impl->e_0)
	{
	}

	ExecutorFacade::~ExecutorFacade()
	{
	}

	void ExecutorFacade::init()
	{
		impl->init();
	}

	void ExecutorFacade::download_data(bool ignore_errors)
	{
		impl->download_data(ignore_errors);
	}

	double ExecutorFacade::time() const
	{
		return impl->time();
	}

	void ExecutorFacade::loop(double* cputimeout, double* gputimeout)
	{
		impl->loop(cputimeout, gputimeout);
	}

	void ExecutorFacade::finish()
	{
		impl->finish();
	}

	void ExecutorFacade::add_job(const std::function<void()>& job)
	{
		impl->add_job(job);
	}
}
}
