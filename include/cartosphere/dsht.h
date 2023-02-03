
#ifndef __DSHT_H__
#define __DSHT_H__

#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

	// This type must be allocated using FFTW's allocators
	// Remove if the following typedef was already taken
	typedef double fftw_real;

	// Discrete spherical harmonic transform
	void cs_fds2ht(size_t B, fftw_real* data, fftw_real* harmonics);

	// Inverse discrete spherical harmonic transform
	void cs_ids2ht(size_t B, fftw_real* harmonics, fftw_real* data);

	// Allocate a workspace for bandlimit B
	double* cs_make_ws2(size_t B);

#ifdef __cplusplus
}
#endif

#endif // !__DSHT_H__
