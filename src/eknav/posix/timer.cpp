#include "eknav/posix/timer.hpp"
#include <cstdlib>

timer::timer()
{
}

void
timer::start()
{
	gettimeofday(&start_t, NULL);
}

double
timer::stop()
{
	struct timeval stop;
	gettimeofday(&stop, NULL);
	if (stop.tv_usec < start_t.tv_usec) {
		stop.tv_usec += 1000000;
		stop.tv_sec -= 1;
	}
	return double(stop.tv_sec - start_t.tv_sec) + (stop.tv_usec - start_t.tv_usec) * 1e-6;
}
