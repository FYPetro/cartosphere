
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

int
cs_index2(int B, int l, int m)
{
	// Format of harmonics: (l,m) with l the degree, m the order
	//
	//     For each row in each triangular part, m is constant,
	//     and as we move right in each row, l increments by 1.
	//     |<-------------- B items, row-major -------------->|
	//   - +--------------------------------------------------+
	//   ^ |_h(__0,_____0)_  h(1, 0)      ...      h(B-1,  0) | UPPER
	//   | | h(B-1,-(B-1)) |_h(1,_1)_     ...      h(B-1,  1) | TRIANGLE
	//   B |      ...          ...   |____...____      ...    | IS FOR
	//   v | h(  1,    -1)     ...     h(B-1,-1) | h(B-1,B-1) | m >= 0
	//   - +-------------------------------------+------------+
	//       STRICTLY LOWER TRIANGLE IS FOR m < 0

	// The upper triangular part is indexed like a regular square
	if (m >= 0)
	{
		return B * m + l;
	}
	// The strictly lower triangular part is indexed like a parallelogram
	// For the row m: there are (B-|m|) rows above it, and degree begins at |m|
	else
	{
		m = -m;
		return B * (B - m) + (l - m);
	}
}

void
cs_fds2ht(int B, double* data, double* harmonics, double* ws2)
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
cs_ids2ht(int B, double* harmonics, double* data, double* ws2,
	fftw_real* scratchpad, fftw_plan idct, fftw_plan idst)
{
	int N = 2 * B;

	// if (B == 2)
	// {
	// 	std::cout << "\nDEBUG: print harmonics per format\n";
	// 	for (int m = 0; m < B; ++m)
	// 	{
	// 		for (int l = B - m; l < B; ++l)
	// 		{
	// 			std::cout << "  "
	// 				<< "B_{" << l << "," << (m - B) << "} = "
	// 				<< harmonics[cs_index2(B, l, (m - B))];
	// 		}
	// 		for (int l = m; l < B; ++l)
	// 		{
	// 			std::cout << "  "
	// 				<< "B_{" << l << "," << m << "} = "
	// 				<< harmonics[cs_index2(B, l, m)];
	// 		}
	// 		std::cout << "\n";
	// 	}
	// }

	// Clear output data and the entire scratchpad
	memset(data, 0, N * N * sizeof(double));
	memset(scratchpad, 0, N * N * 2 * sizeof(double));

	// Compute 1D fourier coefficients for the northern hemisphere first
	// Cosine and sine coefficients are interwoven in the same matrix!
	for (int j = 0; j < N; ++j)
	{
		fftw_real* amj = scratchpad + (2 * N * j);
		fftw_real* bmj = amj + N;
		// Retrieve renormalized P_{l,m} per x_{j}-file
		// This file is already in upper triangular form
		auto rePlmCos = cs_ws2_rePlmCosFile(B, j, ws2);
		// Compute the cosine coefficients
		for (int m = 0; m < B; ++m, ++amj)
		{
			// Compute element-wise product between...
			// 1: ROW m of UPPER TRIANGLE of HARMONICS
			// 2: ROW m of UPPER TRIANGLE rePlmCosFile for x_{j}
			auto row1 = harmonics + (B * m);
			auto row2 = rePlmCos + (2 * B - m + 1) * m / 2;
			for (int l = m; l < B; ++l)
			{
				*amj += row1[l] * row2[l - m];
			}
		}
		// Compute the sine coefficients
		for (int m = 1; m < B; ++m, ++bmj)
		{
			// Compute element-wise product between...
			// 1: ROW B-m of LOWER TRIANGLE of HARMONICS, shifted by m
			// 2: ROW   m of UPPER TRIANGLE rePlmCosFile for x_{j}
			auto row1 = harmonics + (B * (B - m)) - m;
			auto row2 = rePlmCos + (2 * B - m + 1) * m / 2;
			for (int l = m; l < B; ++l)
			{
				*bmj += row1[l] * row2[l - m];
			}
		}
		// Zero out the final sine coefficient
		*bmj++ = 0;
	}

	// if (B == 2)
	// {
	// 	std::cout << "\nDEBUG: DCT-III coefficients\n";
	// 	for (int j = 0; j < N; ++j)
	// 	{
	// 		for (int m = 0; m < B; ++m)
	// 		{
	// 			std::cout << "  "
	// 				<< "a_{" << j << "," << m << "} = "
	// 				<< scratchpad[2 * N * j + m];
	// 		}
	// 		std::cout << "\n";
	// 	}
	// }
	// if (B == 2)
	// {
	// 	std::cout << "\nDEBUG: DST-III coefficients\n";
	// 	for (int j = 0; j < N; ++j)
	// 	{
	// 		for (int m = 0; m < B; ++m)
	// 		{
	// 			std::cout << "  "
	// 				<< "b_{" << j << "," << m << "} = "
	// 				<< scratchpad[2 * N * j + N + m];
	// 		}
	// 		std::cout << "\n";
	// 	}
	// 	std::cout << "\n";
	// }

	// Account for normalization (FFTW to C++17)
	for (int j = 0; j < N; ++j)
	{
		// All DCT-III coefficients with m != 0 must be divided by 2
		auto* target = scratchpad + (2 * N * j + 1);
		for (int k = 1; k < B; ++k)
		{
			*target++ *= 0.5;
		}
		// All DCT-III coefficients with m != B-1 must be divided by 2
		// Note that the final coefficient is zero, so it need not be modified
		target += B;
		for (int k = 0; k < B - 1; ++k)
		{
			*target++ *= 0.5;
		}
	}
	// Perform D{C,S}T-III
	fftw_execute(idct); fftw_execute(idst);
	// Copy results to the eastern hemisphere
	for (int j = 0; j < N; ++j)
	{
		// Aggregate data due to DCT-III
		auto* target = data + (N * j);
		auto* source = scratchpad + (2 * N * j + B);
		for (int k = 0; k < B; ++k)
		{
			*target++ += *source++;
		}
		// Aggregate data due to DST-III
		target = data + (N * j);
		source += B;
		for (int k = 0; k < B; ++k)
		{
			*target++ += *source++;
		}
	}

	// if (B == 2)
	// {
	// 	std::cout << "\nDEBUG: cosine contributions\n";
	// 	for (int j = 0; j < N; ++j)
	// 	{
	// 		for (int m = 0; m < B; ++m)
	// 		{
	// 			std::cout << "  "
	// 				<< "A_{" << j << "," << m << "} = "
	// 				<< scratchpad[2 * N * j + B + m];
	// 		}
	// 		std::cout << "\n";
	// 	}
	// }
	// if (B == 2)
	// {
	// 	std::cout << "\nDEBUG: sine contributions\n";
	// 	for (int j = 0; j < N; ++j)
	// 	{
	// 		for (int m = 0; m < B; ++m)
	// 		{
	// 			std::cout << "  "
	// 				<< "B_{" << j << "," << m << "} = "
	// 				<< scratchpad[2 * N * j + N + B + m];
	// 		}
	// 		std::cout << "\n";
	// 	}
	// 	std::cout << "\n";
	// }

	// Tune coefficients for the western hemisphere
	for (int j = 0; j < N; ++j)
	{
		// For every two DCT-III columns, negate the second
		auto* target = scratchpad + (2 * N * j + 1);
		for (int k = 0; k < B; k += 2, target += 2)
		{
			*target *= -1.0;
		}
		// For every two DST-III columns, negate the first
		target += B - 1;
		for (int k = B; k < N; k += 2, target += 2)
		{
			*target *= -1.0;
		}
	}
	// Perform D{C,S}T-III
	fftw_execute(idct); fftw_execute(idst);
	// Copy results to the western hemisphere
	for (int j = 0; j < N; ++j)
	{
		// Aggregate data due to DCT-III
		auto* target = data + (N * j + B);
		auto* source = scratchpad + (2 * N * j + B);
		for (int k = B; k < N; ++k)
		{
			*target++ += *source++;
		}
		// Aggregate data due to DST-III
		target = data + (N * j + B);
		source += B;
		for (int k = B; k < N; ++k)
		{
			*target++ += *source++;
		}
	}

	// if (B == 2)
	// {
	// 	std::cout << "\nDEBUG: synthesized data\n";
	// 	for (int j = 0; j < N; ++j)
	// 	{
	// 		for (int k = 0; k < N; ++k)
	// 		{
	// 			std::cout << "  u_{" << j << "," << k << "} = " << data[N * j + k];
	// 		}
	// 		std::cout << "\n";
	// 	}
	// 	std::cout << "\n";
	// }
}

void
cs_make_plans2(int B,
	fftw_real* scratchpad, fftw_plan* ptr_idct, fftw_plan* ptr_idst)
{
	int N = 2 * B;

	// For each D{C,S}T-III...
	// Perform rank-1 (1-dimensional)
	int rank = 1;
	// ... of input length B
	int n[] = { B };
	// ... for N batches
	int howmany = { N };

	// The first input element is at
	fftw_real* in = scratchpad;
	// The input of each batch is n-shaped (hence NULL)
	int* inembed = NULL;
	// The stride of each element within each input batch
	int istride = 1;
	// The stride of the first input element across all batches
	int idist = 2 * N;

	// The first output element is at
	fftw_real* out = in + B;
	// The output of each batch is n-shaped (hence NULL)
	int* onembed = NULL;
	// The stride of each element within each output batch
	int ostride = 1;
	// The stride of the first output element across all batches
	int odist = 2 * N;

	// Perform DCT-III for each batch
	fftw_r2r_kind kind[] = { FFTW_REDFT01 };
	// Default runtime flags
	auto flags = FFTW_ESTIMATE;

	// Create DCT-III plan using the advanced real-to-real interface
	*ptr_idct = fftw_plan_many_r2r(rank, n, howmany,
		in, inembed, istride, idist,
		out, onembed, ostride, odist,
		kind, flags);

	// Create DST-III plan using the advanced real-to-real interface
	in += N; out += N; kind[0] = FFTW_RODFT01;
	*ptr_idst = fftw_plan_many_r2r(rank, n, howmany,
		in, inembed, istride, idist,
		out, onembed, ostride, odist,
		kind, flags);
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

	// if (B == 2)
	// {
	// 	std::cout << "\nDEBUG: cs_make_ws2 workspace block 5\n";
	// 	for (int l = 0; l < B; ++l)
	// 	{
	// 		for (int m = 0; m <= l; ++m)
	// 		{
	// 			double* Plms = cs_ws2_rePlmCosRank(B, l, m, ws2);
	// 			std::cout << "~P_{" << l << "," << m << "} =\n";
	// 			for (int j = 0; j < N; ++j)
	// 			{
	// 				std::cout << "  " << Plms[j];
	// 			}
	// 			std::cout << "\n";
	// 		}
	// 	}
	// }

	// [Block 6] Transposed table
	// Memory intensive (lots of strides), consider generating this from scratch
	double** samples = (double **)malloc(B * (B + 1) / 2 * sizeof(double *));
	{
		// Generate N-striding pointers
		auto sourcePtr = blocks[5];
		auto sample = samples;
		for (int m = 0; m < B; ++m)
		{
			for (int l = m; l < B; ++l, sourcePtr += N)
			{
				*sample++ = sourcePtr;
			}
		}
		// For each j, fetch from N-striding pointers
#pragma omp parallel for if (B >= 128)
		for (int j = 0; j < N; ++j)
		{
			auto target = cs_ws2_rePlmCosFile(B, j, ws2);
			auto fileSample = samples;
			for (int m = 0; m < B; ++m)
			{
				for (int l = m; l < B; ++l)
				{
					*(target++) = **(fileSample++);
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

	// if (B == 2)
	// {
	// 	std::cout << "\nDEBUG: cs_make_ws2 workspace block 6\n";
	// 	double* Plms = cs_ws2_rePlmCosFile(B, 0, ws2);
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

	double* file = cs_ws2_rePlmCosRank(B, B, 0, ws2) + (B * (B + 1) / 2 * j);

	return file;
}
