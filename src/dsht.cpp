
#include "cartosphere/dsht.hpp"

#include <xmmintrin.h>

#include <omp.h>

#include "cartosphere/utility.hpp"

// MacOS: pull legendre from boost::math
#if defined(APPLE_LIKE)
#include <boost/math/special_functions/legendre.hpp>
constexpr auto legendre = [](auto &&...args) {
	return boost::math::legendre_p(std::forward<decltype(args)>(args)...);
};
#else
using std::legendre;
#endif

void
cs_fds2ht(int B, fftw_real* data, fftw_real* harmonics, double* ws2)
{
	int N = 2 * B;

	// Perform DCT to resolve hats for m>=0
	// {
	// 	int rank = 2;
	// 	int n[] = { N, 1 };
	// 	int howmany = N;
	//
	// 	fftw_real* in;
	// 	int* inembed = NULL;
	// 	int istride = 1;
	// 	int idist = N;
	//
	// 	fftw_real* out;
	// 	int* onembed = NULL;
	// 	int ostride = 1;
	// 	int odist = N;
	//
	// 	fftw_r2r_kind kind[] = { FFTW_REDFT11 };
	// 	auto flags = FFTW_MEASURE;
	// 	fftw_plan_many_r2r(rank, n, howmany, in, inembed, istride, idist,
	// 		out, onembed, ostride, odist, kind, flags);
	// }
	// // Perform DST to resolve hats for m<0
}

void
cs_ids2ht(int B, double* harmonics, double* data, double* ws2)
{
	int N = 2 * B;
}

double*
cs_make_ws2(int B)
{
	int N = 2 * B;

	// Allocate workspace
	double* const ws2 = new double
		[(4 + 3 * N + (N - 1) * N + N * B * (B + 1))];

	// Segmentize the workspace into multiple blocks
	double* const blocks[7] = {
		// Block 0: 4 elements
		// Element 0: bandlimit
		// Element 1-3: unused
		ws2,

		// Block 1: N elements
		// Stores the weights for each cell through w = P\(2pi/B delta)
		ws2 + 4,

		// Block 2: N elements
		// Stores the polar cosines: cos(theta_{j}) = x_{j}
		ws2 + (4 + N),

		// Block 3: N elements
		// Stores the polar sines: sin(theta_{j}) = y_{j} = sqrt(1-x^2)
		ws2 + (4 + 2 * N),

		// Block 4: (N-1)*N elements
		// Stores the azimuthal sines and cosines for each azimuth
		ws2 + (4 + 3 * N), // sin(m phi_{k}) ... cos(m phi_{k}) for each k

		// Block 5: N*B*(B+1)/2 elements
		// Stores the C++17 renormalized ~P_{l,m} = q_{l}^{m}P_{l}^{m}(x_{j})
		// Dimensions: First j, then m, then l
		ws2 + (4 + 3 * N + (N - 1) * N),

		// Block 6: N*B*(B+1)/2 elements
		// Permutes the block above
		ws2 + (4 + 3 * N + (N - 1) * N + N * B * (B + 1) / 2)
	};
	
	// [Block 0] Bandlimit
	blocks[0][0] = B;

	// [Block 1-2] Generate weights by solving, for 0 <= l < N = 2B,
	// \sum_{j=1}^{N}(P_{l}(cos(theta_{j})))w_{j}=(2pi/B)delta_{0,l}
	double* w = blocks[1];
	double* x = blocks[2];
	double* tempCosPls = blocks[3]; // Do not overwrite until moved!
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
				target[j] = legendre(l, x[j]);
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

	// [Block 3, 5] Populate associated Legendre table recursively
	double* y = blocks[3];
	double* reCosPlms = blocks[5];
	{
		// Move tempCosPls to correct location.
		for (int l = B - 1; l >= 0; --l)
		{
			double* source = tempCosPls + (N * l);
			double* target = cs_ws2_rePlmCosRank(B, l, 0, ws2);
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
		double* Plls = cs_ws2_rePlmCosRank(B, 0, 0, ws2);
		for (int l = 0; l < B - 1; ++l)
		{
			double* Plp1_lp1s = cs_ws2_rePlmCosRank(B, l + 1, l + 1, ws2);
			for (int j = 0; j < N; ++j)
			{
				Plp1_lp1s[j] = (2 * l + 1) * y[j] * Plls[j];
			}
			double* offDiagonal = cs_ws2_rePlmCosRank(B, l + 1, l + 1, ws2);
			Plls = Plp1_lp1s;
		}

		// Populate off-diagonal: P_{l}^{l} => P_{l+1}^{l}
		for (int l = 1; l < B - 1; ++l)
		{
			Plls = cs_ws2_rePlmCosRank(B, l, l, ws2);

			double* Plp1_ls = cs_ws2_rePlmCosRank(B, l + 1, l, ws2);
			for (int j = 0; j < N; ++j)
			{
				Plp1_ls[j] = (2 * l + 1) * x[j] * Plls[j];
			}
		}

		// Populate horizontally: P_{l-1}^{m} & P_{l}^{m} => P_{l+1}^{m}
#pragma omp parallel for if (B >= 128)
		for (int m = 1; m < B - 1; ++m)
		{
			double* Plm1_ms = cs_ws2_rePlmCosRank(B, m, m, ws2);
			double* Plms = cs_ws2_rePlmCosRank(B, m + 1, m, ws2);
			for (int l = m + 1; l < B - 1; ++l)
			{
				double* Plp1_ms = cs_ws2_rePlmCosRank(B, l + 1, m, ws2);
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

				double* Plms = cs_ws2_rePlmCosRank(B, l, m, ws2);
				for (int j = 0; j < N; ++j)
				{
					Plms[j] *= qlm;
				}
			}
		}
	}

	// [Block 4] Populate trig values for inverse transform
	double* trigs = blocks[4];
	{
		// Cosine and sine values of the same azimuth are consecutive
		double* ptr = trigs;

		// Loop through each azimuth
		for (int k = 0; k < N; ++k)
		{
			double phi = 2 * M_PI * ((k + 0.5) / N);
			
			// -B < m < 0
			for (int m = B - 1; m > 0; --m)
			{
				*ptr++ = sin(m * phi);
			}

			// m = 0
			*ptr++ = 0;

			// 0 < m < B
			for (int m = 1; m < B; ++m)
			{
				*ptr++ = cos(m * phi);
			}
		}
	}

	// Debug Block 5
	// if (B == 4)
	// {
	// 	std::cout << "\nDEBUG:\n";
	// 	for (int l = 0; l < B; ++l)
	// 	{
	// 		for (int m = 0; m <= l; ++m)
	// 		{
	// 			double* Plms = cs_ws2_rePlmRank(B, l, m, ws2);
	//
	// 			std::cout << "~P_{" << l << "," << m << "} =\n";
	// 			for (int j = 0; j < N; ++j)
	// 			{
	// 				std::cout << "  " << Plms[j];
	// 			}
	// 			std::cout << "\n";
	// 		}
	// 	}
	// 	std::cout << "\n";
	// }

	// [Block 6] Transposed table
	// Memory intensive (lots of strides), consider generating this from scratch
	double** samples = (double **)malloc(B * (B + 1) / 2 * sizeof(double *));
	{
		// Generate N-striding pointers
		auto sourcePtr = blocks[5];
		auto stridingPtr = samples;
		for (int m = 0; m < B; ++m)
		{
			for (int l = m; l < B; ++l)
			{
				*(stridingPtr++) = sourcePtr;
				sourcePtr += N;
			}
		}
		// For each j, fetch from N-striding pointers
#pragma omp parallel for if (B >= 128)
		for (int j = 0; j < N; ++j)
		{
			auto target = cs_ws2_rePlmCosFile(B, j, ws2);
			auto stridingPtrPtr = samples;
			for (int m = 0; m < B; ++m)
			{
				for (int l = m; l < B; ++l)
				{
					*(target++) = **(stridingPtrPtr++);
				}
			}
			// Shift all N-striding pointers by 1
			for (int i = 0; i < B * (B + 1) / 2; ++i)
			{
				samples[i] += 1;
			}
		}
	}
	free(samples);
	samples = nullptr;

	// Debug Block 6
	// if (B == 4)
	// {
	// 	std::cout << "\nDEBUG:\n";
	// 	double* Plms = cs_ws2_rePlmFile(B, 0, ws2);
	// 	for (int j = 0; j < N; ++j)
	// 	{
	// 		std::cout << "For phi_{" << j << "}:\n";
	// 		for (int m = 0; m < B; ++m)
	// 		{
	// 			for (int l = m; l < B; ++l)
	// 			{
	// 				std::cout << "  ~P_{" << l << "," << m << "} = " << *Plms++;
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
	free(ws2);
}

double*
cs_ws2_rePlmCosRank(int B, int l, int m, double* ws2)
{
	int N = 2 * B;

	double* rank = ws2 +
		(4 + 3 * N + (N - 1) * N + N * (l * (l + 1) / 2 + m));

	return rank;
}

double*
cs_ws2_rePlmCosFile(int B, int j, double* ws2)
{
	int N = 2 * B;

	double* file = cs_ws2_rePlmCosRank(B, B, 0, ws2) + j;

	return file;
}
