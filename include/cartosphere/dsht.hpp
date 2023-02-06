
#ifndef __DSHT_H__
#define __DSHT_H__

#include <fftw3.h>

// This type must be allocated using FFTW's allocators
// Remove if the following typedef was already taken
typedef double fftw_real;

// Discrete spherical harmonic transform
void cs_fds2ht(int B, double* data, double* harmonics, double* ws2);

// Inverse discrete spherical harmonic transform
void cs_ids2ht(int B, double* harmonics, double* data, double* ws2);

// Allocate a workspace for bandlimit B
double* cs_make_ws2(int B);

// Destroy workspace
void cs_free_ws2(double* ws2);

// Fetch
double* cs_ws2_rePlmCosRank(int B, int l, int m, double* ws2);

// Fetch
double* cs_ws2_rePlmCosFile(int B, int j, double* ws2);

#endif // !__DSHT_H__
