#include "random_seed.hpp"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <wincrypt.h>

uint32_t random_seed()
{
	HCRYPTPROV provider;
	CryptAcquireContext(&provider, NULL, NULL, PROV_RSA_AES, 0);

	uint32_t ret;
	CryptGenRandom(provider, sizeof(ret), (BYTE*)&ret);
	CryptReleaseContext(provider, 0);
	return ret;
}

