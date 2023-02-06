
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
	double* const ws2 = new double[(1 + 3 * N + N * B * (B + 1) / 2)];

	double* const blocks[5] = {
		ws2, // B is the bandlimit
		ws2 + 1, // w = P_{l}^{m}(x)\[2pi/B; 0; ... 0]
		ws2 + (1 + N), // cos(theta_{j}) = x_{j}
		ws2 + (1 + 2 * N), // sin(theta_{j}) = y_{j} = sqrt(1-x^2)
		ws2 + (1 + 3 * N), // q_{l}^{m}P_{l}^{m}(x_{j})
	};
	
	// [Block 0] Bandlimit
	blocks[0][0] = B;

	// [Block 1-2] Generate weights by solving, for 0 <= l < N = 2B,
	// \sum_{j=1}^{N}(P_{l}(cos(theta_{j})))w_{j}=(2pi/B)delta_{0,l}
	double* w = blocks[1];
	double* x = blocks[2];
	double* tempCosPls = blocks[3];
	{
		// Compute the cosine of polar angles
		for (int j = 0; j < N; ++j)
		{
			x[j] = cos(M_PI / N * (j + 0.5));
		}

		// Compute weights form Legendre coefficients up to degree 2B-1
#pragma omp parallel for if (B >= 128)
		for (int l = 0; l < N; ++l)
		{
			double* target = tempCosPls + (N * l);
			for (int j = 0; j < N; ++j)
			{
				target[j] = std::legendre(l, x[j]);
			}
		}

		// Fill Eigen matrices A, b, solve Au=b, extract results
		Matrix A = Eigen::Map<Matrix>(tempCosPls, N, N);
		Vector b(N);
		memset(b.data(), 0, N * sizeof(double));
		b[0] = 2 * M_PI / B;
		
		Vector u = A.partialPivLu().solve(b); // PartialPivLU suffices
		memcpy(w, u.data(), N * sizeof(double));
	}
	
	// [Block 3-4] Populate associated Legendre table recursively
	double* y = blocks[3];
	double* reCosPlms = blocks[4];
	{
		// Move tempCosPls to correct location.
		for (int l = B - 1; l >= 0; --l)
		{
			double* source = tempCosPls + (N * l);
			double* target = cs_ws2_reCosPlms(B, l, 0, ws2);
			for (int j = 0; j < N; ++j)
			{
				memcpy(target, source, N * sizeof(double));
			}
		}
		tempCosPls = nullptr;

		// Compute the abs(sin(@)) of polar angles
		for (int j = 0; j < N; ++j)
		{
			y[j] = sin(M_PI / N * (j + 0.5));
		}

		// Populate diagonal: P_{l}^{l} => P_{l+1}^{l^1}
		double* Plls = cs_ws2_reCosPlms(B, 0, 0, ws2);
		for (int l = 0; l < B - 1; ++l)
		{
			double* Plp1_lp1s = cs_ws2_reCosPlms(B, l + 1, l + 1, ws2);
			for (int j = 0; j < N; ++j)
			{
				Plp1_lp1s[j] = (2 * l + 1) * y[j] * Plls[j];
			}
			double* offDiagonal = cs_ws2_reCosPlms(B, l + 1, l + 1, ws2);
			Plls = Plp1_lp1s;
		}

		// Populate off-diagonal: P_{l}^{l} => P_{l+1}^{l}
		for (int l = 1; l < B - 1; ++l)
		{
			Plls = cs_ws2_reCosPlms(B, l, l, ws2);

			double* Plp1_ls = cs_ws2_reCosPlms(B, l + 1, l, ws2);
			for (int j = 0; j < N; ++j)
			{
				Plp1_ls[j] = (2 * l + 1) * x[j] * Plls[j];
			}
		}

		// Populate horizontally: P_{l-1}^{m} & P_{l}^{m} => P_{l+1}^{m}
#pragma omp parallel for if (B >= 128)
		for (int m = 1; m < B - 1; ++m)
		{
			double* Plm1_ms = cs_ws2_reCosPlms(B, m, m, ws2);
			double* Plms = cs_ws2_reCosPlms(B, m + 1, m, ws2);
			for (int l = m + 1; l < B - 1; ++l)
			{
				double* Plp1_ms = cs_ws2_reCosPlms(B, l + 1, m, ws2);
				for (int j = 0; j < N; ++j)
				{
					Plp1_ms[j] = ((2 * l + 1) * x[j] * Plms[j] - (l + m) * Plm1_ms[j])
						/ (l - m + 1);
				}
				Plm1_ms = Plms;
				Plms = Plp1_ms;
			}
		}

		// Divide stuff
#pragma omp parallel for if (B >= 128)
		for (int l = 0; l < B; ++l)
		{
			for (int m = 0; m <= l; ++m)
			{
				double qlm = sqrt((l + 0.5) / M_PI);
				for (int i = l - m + 1; i <= l + m; ++i)
				{
					qlm /= sqrt(i);
				}

				double* Plms = cs_ws2_reCosPlms(B, l, m, ws2);
				for (int j = 0; j < N; ++j)
				{
					Plms[j] *= qlm;
				}
			}
		}
	}

	// if (B == 4)
	// {
	// 	std::cout << "\nDEBUG:\n";
	// 	for (int l = 0; l < B; ++l)
	// 	{
	// 		for (int m = 0; m <= l; ++m)
	// 		{
	// 			double* Plms = cs_ws2_reCosPlms(B, l, m, ws2);
	//
	// 			std::cout << "~P_{" << l << "}_{" << m << "} =\n";
	// 			for (int j = 0; j < N; ++j)
	// 			{
	// 				std::cout << "  " << Plms[j];
	// 			}
	// 			std::cout << "\n";
	// 		}
	// 	}
	// 	std::cout << "\n";
	// }

	return ws2;
}

void
cs_free_ws2(double *ws2)
{
	delete[] ws2;
}

double*
cs_ws2_reCosPlms(int B, int l, int m, double* ws2)
{
	int N = 2 * B;

	double* reCosPlms = ws2 + (1 + 3 * N + N * (l * (l + 1) / 2 + m));

	return reCosPlms;
}
