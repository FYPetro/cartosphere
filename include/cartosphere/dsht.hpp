
#ifndef __DSHT_H__
#define __DSHT_H__

#include <fftw3.h>

// This type must be allocated using FFTW's allocators
// Remove if the following typedef was already taken
typedef double fftw_real;

// Generate the linear index for degree l, order m in bandlimit-B harmonics
int cs_index2(int B, int l, int m);

// Discrete spherical harmonic transform
void cs_fds2ht(int B, double* data, double* harmonics, double* ws2);

// Inverse discrete spherical harmonic transform
void cs_ids2ht(int B, double* harmonics, double* data, double* ws2,
	fftw_real* scratchpad, fftw_plan dct3, fftw_plan dst3);

// Generate, semi-interweaved DCT-III and DST-III plans for cs_ids2ht usage
//      // Assume harmonics is B * B and data is N * N
//      int N = 2 * B;
//      fftw_real* scratchpad = fftw_alloc_real(N * N * 2);
//      // Create Type III plans
//      fftw_plan idct, idst;
//      cs_make_plans2(B, scratchpad, &idct, &idst);
//      cs_ids2ht(B, harmonics, data, ws2, scratchpad, idct, idst);
//      // Deallocate plans and scratchpad when no longer needed
//      fftw_destroy_plan(idct);
//      fftw_destroy_plan(idst);
//      fftw_free(scratchpad)
void cs_make_plans2(int B,
	fftw_real* scratchpad, fftw_plan *ptr_idct, fftw_plan *ptr_idst);

// Allocate a workspace for bandlimit B
// WARNING: B must be a positive even number!
// WARNING: For best performance, B must be a power of 2!
// No current plan to work with odd bandlimits, because that's just odd!
double* cs_make_ws2(int B);

// Destroy workspace
void cs_free_ws2(double* ws2);

// Fetch
double* cs_ws2_rePlmCosRank(int B, int l, int m, double* ws2);

// Fetch
double* cs_ws2_rePlmCosFile(int B, int j, double* ws2);

#endif // !__DSHT_H__
