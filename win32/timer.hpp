#ifndef AHRS_TIMER_HPP
#define AHRS_TIMER_HPP

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

class timer
{
	static LARGE_INTEGER clockspeed;
	LARGE_INTEGER last_start;
	LARGE_INTEGER last_stop;

public:
	timer();
	void start();
	double stop();
};
#endif
