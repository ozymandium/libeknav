#include "timer.hpp"

LARGE_INTEGER timer::clockspeed;

timer::timer()
{
	if (clockspeed.QuadPart == 0) {
		QueryPerformanceFrequency(&clockspeed);
	}
	last_start.QuadPart = 0;
}

void timer::start()
{
	QueryPerformanceCounter(&last_start);
}

double timer::stop()
{
	LARGE_INTEGER stop;
	QueryPerformanceCounter(&stop);
	return double(stop.QuadPart - last_start.QuadPart) / clockspeed.QuadPart;
}
