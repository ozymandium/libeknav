#include "posix/random_seed.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

uint32_t random_seed()
{
	int fd = open("/dev/urandom", O_RDONLY);
	uint32_t ret;
	read(fd, &ret, sizeof(ret));
	close(fd);
	return ret;
}
