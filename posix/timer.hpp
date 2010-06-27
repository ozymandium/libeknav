#ifndef AHRS_TIMER_HPP
#define AHRS_TIMER_HPP

#include <sys/time.h>

class timer
{
	struct timeval start_t;
public:
	timer();
	void start();
	double stop();
};
#endif
