
#ifndef __DSHT_H__
#define __DSHT_H__

#include <fftw3.h>

// This type must be allocated using FFTW's allocators
// Remove if the following typedef was already taken
typedef double fftw_real;

// Discrete spherical harmonic transform
void cs_fds2ht(int B, double* data, double* harmonics, double* ws2);

// Inverse discrete spherical harmonic transform
void cs_ids2ht(int B, double* harmonics, double* data, double* ws2,
	fftw_real* scratchpad, fftw_plan plan);

// Generate, from user allocated fftw-ready arrays, an interweaved half DCT-III,
// half DST-III many real-to-real plan. Example snippet below:
//      // Assume harmonics is B * B double and data is N * N
//      int N = 2 * B;
//      fftw_real* scratchpad = fftw_malloc(N * N * 2);
//      auto plan = cs_ids2ht_plan(B, scratchpad); // Create said plan
//      cs_ids2ht(B, harmonics, data, ws2, scratchpad, plan);
//      fftw_destroy_plan(plan); // Deallocate the plan when no longer needed
fftw_plan cs_ids2ht_plan(int B, fftw_real* scratchpad);

// Allocate a workspace for bandlimit B
// WARNING: B must be a positive even numbers!
// WARNING: For best performance, B must be a power of 2!
// There is currently no plan to extend B into all positive integers
double* cs_make_ws2(int B);

// Destroy workspace
void cs_free_ws2(double* ws2);

// Fetch
double* cs_ws2_rePlmCosRank(int B, int l, int m, double* ws2);

// Fetch
double* cs_ws2_rePlmCosFile(int B, int j, double* ws2);

#endif // !__DSHT_H__
