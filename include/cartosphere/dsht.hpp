
#ifndef __DSHT_H__
#define __DSHT_H__

#include <fftw3.h>

// This type must be allocated using FFTW's allocators
// Remove if the following typedef was already taken
typedef double fftw_real;

// Generate the linear index for degree l, order m in bandlimit-B harmonics
int cs_index2(int B, int l, int m);

// Discrete spherical harmonic transform
void cs_fds2ht(int B, double* data, double* harmonics, double* ws2,
	fftw_real* scratchpad, fftw_plan many_dct, fftw_plan many_dst);

// Inverse discrete spherical harmonic transform
void cs_ids2ht(int B, double* harmonics, double* data, double* ws2,
	fftw_real* scratchpad, fftw_plan many_idct, fftw_plan many_idst);

// Generate, semi-interweaved DCT-II and DST-II plans for cs_fds2ht usage
void cs_ids2ht_plans(int B,
	fftw_real* scratchpad, fftw_plan* ptr_many_dct, fftw_plan* ptr_many_dst);

// Generate, semi-interweaved DCT-III and DST-III plans for cs_ids2ht usage
//      // Assume harmonics is B * B and data is N * N
//      int N = 2 * B;
//      fftw_real* scratchpad = fftw_alloc_real(N * N * 2);
//      // Create Type III plans using the fftw_plan_many_r2r interface
//      fftw_plan many_idct, many_idst;
//      cs_ids2ht_plans(B, scratchpad, &many_idct, &many_idst);
//      cs_ids2ht(B, harmonics, data, ws2, scratchpad, many_idct, many_idst);
//      // Deallocate plans and scratchpad when no longer needed
//      fftw_destroy_plan(many_idct);
//      fftw_destroy_plan(many_idst);
//      fftw_free(scratchpad)
void cs_ids2ht_plans(int B,
	fftw_real* scratchpad, fftw_plan* ptr_many_idct, fftw_plan* ptr_many_idst);

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
