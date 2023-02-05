
#include "cartosphere/dsht.hpp"

#include <xmmintrin.h>

#include <omp.h>

#include "cartosphere/utility.hpp"

void
cs_fds2ht(int B, fftw_real* data, fftw_real* harmonics, double* ws2)
{
	int N = 2 * B;
}

void
cs_ids2ht(int B, fftw_real* harmonics, fftw_real* data, double* ws2)
{
	int N = 2 * B;

	for (int l = 0; l < B; ++l)
	{
		// Perform DCT to resolve hats for m>=0

		// Perform DST to resolve hats for m<0
	}
}

double*
cs_make_ws2(int B)
{
	int N = 2 * B;

	// Allocate workspace
	double* ws2 = new double[(1 + N + N * N * (N + 1) / 2)];

	// [Block 1] Bandlimit
	ws2[0] = B;

	// [Block 2] Generate weights by solving, for 0 <= l < N = 2B,
	// \sum_{j=1}^{N}(P_{l}(cos(theta_{j})))w_{j}=(2pi/B)delta_{0,l}
	double* w = ws2 + 1;
	{
		double* tempCosPolars = w;
		// Compute the cosine of polar angles
		for (int j = 0; j < N; ++j)
		{
			tempCosPolars[j] = cos(M_PI / N * (j + 0.5));
		}

		// Generate all Legendre coefficients (m=0) and compute weights
		double* tempCosPls = w + N;
#pragma omp parallel for if (B >= 128)
		for (int l = 0; l < N; ++l)
		{
			double* target = tempCosPls + (N * l);
			for (int j = 0; j < N; ++j)
			{
				target[j] = std::legendre(l, tempCosPolars[j]);
			}
		}

		// Fill Eigen matrices A, b, solve for x, extract results
		Matrix A = Eigen::Map<Matrix>(tempCosPls, N, N);
		Vector b(N);
		memset(b.data(), 0, N * sizeof(double));
		b[0] = 2 * M_PI / B;
		
		Vector x = A.partialPivLu().solve(b); // PartialPivLU is good enough
		memcpy(w, x.data(), N * sizeof(double));
	}
	
	// Populate the associated Legendre table recursively

	return ws2;
}

void
cs_free_ws2(double *ws2)
{
	delete[] ws2;
}

double*
cs_ws2_cosPlms(int B, int l, int m, double* ws2)
{
	int N = 2 * B;

	double* cosPlms = ws2 + (1 + N + N * (l * (l + 1) / 2 + m));

	return cosPlms;
}
