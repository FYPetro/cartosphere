
#ifndef __DSHT_HPP__
#define __DSHT_HPP__

#include <fftw3.h>

// This type must be allocated using FFTW's allocators
// Remove if the following typedef was already taken
typedef double fftw_real;

// Generate the linear index for degree l, order m in bandlimit-B harmonics
int cs_index2(int B, int l, int m);

// Generate the linear index for degree l, order m for associated Legendre
// Note that m>=0, index first m, then l
int cs_index2_assoc(int B, int l, int m);

// Discrete spherical harmonic transform
void cs_fds2ht(int B, const double* data, double* harmonics, const double* ws2);

// Inverse discrete spherical harmonic transform
void cs_ids2ht(int B, const double* harmonics, double* data, const double* ws2,
	fftw_real* pad, fftw_plan many_idct, fftw_plan many_idst);

// Partial derivative w.r.t. theta, inverted form harmonics
void cs_ids2ht_dp(int B, const double* harmonics, double* partials, const double* ws2,
	fftw_real* pad, fftw_plan many_idct, fftw_plan many_idst);

// Partial derivative w.r.t. phi, inverted from harmonics
void cs_ids2ht_da(int B, const double* harmonics, double* partials, const double* ws2,
	fftw_real* pad, fftw_plan many_idct, fftw_plan many_idst);

// Generate, semi-interweaved DCT-III and DST-III plans for cs_ids2ht usage
//      // Assume harmonics is B * B and data is N * N
//      int N = 2 * B;
//      fftw_real* pad = fftw_alloc_real(N * N * 2);
//      // Create Type III plans using the fftw_plan_many_r2r interface
//      fftw_plan many_idct, many_idst;
//      cs_ids2ht_plans(B, pad, &many_idct, &many_idst);
//      cs_ids2ht(B, harmonics, data, ws2, pad, many_idct, many_idst);
//      // Deallocate plans and scratchpad when no longer needed
//      fftw_destroy_plan(many_idct);
//      fftw_destroy_plan(many_idst);
//      fftw_free(scratchpad)
void cs_ids2ht_plans(int B, fftw_real* pad,
	fftw_plan* ptr_many_idct, fftw_plan* ptr_many_idst);

// Given a valid scratch pad
// Properly execute FFTW plans to obtain desired inverse transform
// Internal to cs_ids2ht, cs_ids2ht_dp, cs_ids2ht_da
void cs_ids2ht_execute(int B, fftw_real* pad, fftw_real* data,
	fftw_plan many_idct, fftw_plan many_idst);

// Allocate a workspace for bandlimit B
// Remember to free it using delete[]!
// WARNING: B must be a positive even number!
// WARNING: For best performance, B must be a power of 2!
// No current plan to work with odd bandlimits, because that's just odd!
double* cs_make_ws2(int B);

void cs_make_ws2(int B, double* ws2);

// Returns the size of the workspace
int cs_ws2_size(int B);

// Fetch
double* cs_ws2_rePlmCosRank(int B, int l, int m, double* ws2);
const double* cs_ws2_rePlmCosRank(int B, int l, int m, const double* ws2);

// Fetch
double* cs_ws2_rePlmCosFile(int B, int j, double* ws2);
const double* cs_ws2_rePlmCosFile(int B, int j, const double* ws2);

// Fetch
double* cs_ws2_drePlmCosFile(int B, int j, double* ws2);
const double* cs_ws2_drePlmCosFile(int B, int j, const double* ws2);

#endif // !__DSHT_H__
