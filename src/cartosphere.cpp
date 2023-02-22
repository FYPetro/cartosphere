
#include "cartosphere/cartosphere.hpp"
using Cartosphere::Point;
using Cartosphere::SpectralGlobe;
using Cartosphere::FiniteElementGlobe;

#include "cartosphere/dsht.hpp"

void
SpectralGlobe::initialize_solver()
{
	// B is treated as the NEW bandlimit
	// N is treated as twice the OLD bandlimit
	int n = B * 2;

	// If B==0, deallocate, reset
	// If B!=0 and n==N, reset, initialize
	// If B!=0 and n!=N, deallocate, allocate, reset, initialize
	if (B == 0 || n != N)
	{
		// Deallocate
		cleanup();
	}
	if (n != N)
	{
		N = n;
		// Resize
		init_data.resize(N * N);
		time_data.resize(N * N);
		init_hats.resize(B * B);
		time_hats.resize(B * B);
		time_dp.resize(N * N);
		time_da.resize(N* N);
		time_grad.resize(N * N);
		// Allocate
		if (B > 0)
		{
			ws2.resize(cs_ws2_size(B));
			cs_make_ws2(B, ws2.data());
			ipad = fftw_alloc_real(N * N * 2);
			cs_ids2ht_plans(B, ipad, &idct, &idst);
		}
	}
	
	// Reset
	history.clear();
	time_data_north = time_data_south = 0;
	time_grad_north = time_grad_south = { 0, 0, 0 };
	
	// Initialize
	if (B > 0)
	{
		// Sample initial condition
		double phi, theta;
		for (int j = 0; j < N; ++j)
		{
			theta = M_PI / N * (j + 0.5);
			for (int k = 0; k < N; ++k)
			{
				phi = M_PI / B * (k + 0.5);
				Point P(theta, phi);
				init_data[N * j + k] = initFunction(P);
			}
		}

		// Compute initial Fourier coefficients
		cs_fds2ht(B, init_data.data(), init_hats.data(), ws2.data());

		// Update initial data and gradient
		advance_solver(0, 0);
	}
}

void
SpectralGlobe::cleanup()
{
	if (ipad != nullptr)
	{
		fftw_free(ipad);
		ipad = nullptr;
	}
	if (idct != NULL)
	{
		fftw_destroy_plan(idct);
		idct = NULL;
	}
	if (idst != NULL)
	{
		fftw_destroy_plan(idst);
		idst = NULL;
	}
}

void
SpectralGlobe::advance_solver(double time, double delta)
{
	// Interval: [0, t]
	double t = time + delta;

	double* H = time_hats.data();
	double* D = time_data.data();
	double* P[2] = { time_dp.data(), time_da.data() };
	double* W = ws2.data();
	
	// Compute decayed coefficients
	{
		int l, m, i;
		for (l = 0; l < B; ++l)
		{
			int eigenvalue = -l * (l + 1);
			for (m = -l; m <= l; ++m)
			{
				i = cs_index2(B, l, m);
				H[i] = init_hats[i] * exp(eigenvalue * t);
			}
		}
	}

	// Compute homogenized data
	cs_ids2ht(B, H, D, W, ipad, idct, idst);

	// Compute a velocity field at each grid cell corner
	cs_ids2ht_dp(B, H, P[0], W, ipad, idct, idst);
	cs_ids2ht_da(B, H, P[1], W, ipad, idct, idst);
	{
		int i = 0;
		double theta, phi;
		double cos_theta, sin_theta, cos_phi, sin_phi;
		// Turn dp e_theta + da e_phi into cartesian coordinates
		for (int j = 0; j < N; ++j)
		{
			theta = M_PI / N * (j + 0.5);
			cos_theta = cos(theta);
			sin_theta = sin(theta);
			for (int k = 0; k < N; ++k, ++i)
			{
				phi = M_PI / B * (k + 0.5);
				cos_phi = cos(phi);
				sin_phi = sin(phi);
				// Convert gradient!
				auto& grad = time_grad[i];
				grad.x = P[0][i] * cos_theta * cos_phi - P[1][i] * sin_phi / sin_theta;
				grad.y = P[0][i] * cos_theta * sin_phi + P[1][i] * cos_phi / sin_theta;
				grad.z = P[0][i] * (-sin_theta);
			}
		}
	}

	// Compute data and velocities at the poles
	time_data_north = 0;
	time_data_south = 0;
	for (int l = 0; l < B; ++l)
	{
		double q_l = sqrt((2 * l + 1) / (4 * M_PI));
		double q_hat_l = q_l * time_hats[cs_index2(B, l, 0)];
		time_data_north += q_hat_l;
		time_data_south += ((l % 2) ? -1 : 1) * q_hat_l;
	}
	time_grad_north = { 0, 0, 0 };
	time_grad_south = { 0, 0, 0 };
	{
		double phi, cos_phi, sin_phi;
		double sin_theta = sin(M_PI / N * 0.5);
		int offset = N * (N - 1);
		for (int k = 0; k < N; ++k)
		{
			phi = M_PI / B * (k + 0.5);
			cos_phi = cos(phi);
			sin_phi = sin(phi);
			time_grad_north.x += P[0][k] * cos_phi - P[1][k] * sin_phi / sin_theta;
			time_grad_north.y += P[0][k] * sin_phi + P[1][k] * cos_phi / sin_theta;
			time_grad_south.x += -P[0][offset + k] * cos_phi - P[1][offset + k] * sin_phi / sin_theta;
			time_grad_south.y += -P[0][offset + k] * sin_phi + P[1][offset + k] * cos_phi / sin_theta;
		}
	}
	time_grad_north /= N;
	time_grad_south /= N;
}

void
SpectralGlobe::velocity(const vector<Point>& points, vector<FL3>& velocities) const
{
	// Prepare for logging
	Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

	if (FLAGS_minloglevel == 0)
	{
		LOG(INFO) << "SpectralGlobe::velocity partial_theta\n"
			<< Eigen::Map<const MatrixRowMajor>(time_dp.data(), N, N).format(OctaveFmt);
		LOG(INFO) << "SpectralGlobe::velocity partial_phi\n"
			<< Eigen::Map<const MatrixRowMajor>(time_da.data(), N, N).format(OctaveFmt);
	}

	if (FLAGS_minloglevel == 0)
	{
		stringstream sst;
		sst << "SpectralGlobe::velocity time_grad\n"
			<< "  NORTH: [" << time_grad_north.x << ", "
			<< time_grad_north.y << ", " << time_grad_north.z << "]\n"
			<< "  SOUTH: [" << time_grad_south.x << ", "
			<< time_grad_south.y << ", " << time_grad_south.z << "]\n";
		for (int j = 0; j < N; ++j)
		{
			sst << "  ";
			for (int k = 0; k < N; ++k)
			{
				auto& g = time_grad[j * N + k];
				sst << "(" << j << "," << k << ") = "
					<< "[" << g.x << ", " << g.y << ", " << g.z << "] ";
			}
			sst << "\n";
		}
		LOG(INFO) << sst.str();
	}

	// Calculate the j, k index of the cell that contains each point
	// Then compute the velocity based on the shape of the cell
	for (int i = 0; i < points.size(); ++i)
	{
		const Point& P = points[i];
		
		// Compute the fractional j, k indices aligned with cell centers
		double j_frac = P.p() * N * M_1_PI - 0.5;
		double k_frac = P.a() * B * M_1_PI - 0.5;
		if (k_frac < 0)
		{
			k_frac += N;
		}
		
		// Compute whole j, k indices aligned with cell centers
		int j_n = (int)floor(j_frac);
		int j_s = (j_n + 1);
		int k_w = (int)floor(k_frac) % N;
		int k_e = (k_w + 1) % N;
		
		// Compute remainder coordinates for later bilinear interpolation
		j_frac -= j_n;
		k_frac -= k_w;

		// Obtain D and gradient on the northern side
		double data_n;
		FL3 grad_n;
		if (j_n == -1)
		{
			// The northern edge is degenerate; it is the north pole
			data_n = time_data_north;
			grad_n = time_grad_north;
			// Scale j_frac from [0.5,1] to [0,1]
			j_frac = 2 * j_frac - 1;
		}
		else
		{
			// The northern edge is part of a latitude circle
			data_n = (1 - k_frac) * time_data[N * j_n + k_w]
				+ k_frac * time_data[N * j_n + k_e];
			grad_n = (1 - k_frac) * time_grad[N * j_n + k_w]
				+ k_frac * time_grad[N * j_n + k_e];
		}

		// Obtain data and gradient on the southern side
		double data_s;
		FL3 grad_s;
		if (j_s == N)
		{
			// The southern edge is degenerate; it is the north pole
			data_s = time_data_south;
			grad_s = time_grad_south;
			// Scale j_frac from [0,0.5] to [0,1]
			j_frac = 2 * j_frac;
		}
		else
		{
			// The southern edge is part of a latitude circle
			data_s = (1 - k_frac) * time_data[N * j_s + k_w]
				+ k_frac * time_data[N * j_s + k_e];
			grad_s = (1 - k_frac) * time_grad[N * j_s + k_w]
				+ k_frac * time_grad[N * j_s + k_e];
		}

		// Perform bilinear interpolation
		double data = (1 - j_frac) * data_n + j_frac * data_s;
		FL3 grad = (1 - j_frac) * grad_n + j_frac * grad_s;

		// Compute the velocity
		velocities[i] = -grad / data;
	}
}

void
FiniteElementGlobe::initialize_solver()
{
	// TODO: Implement
}

void
FiniteElementGlobe::advance_solver(double time, double delta)
{
	// TODO: Implement
}

void
FiniteElementGlobe::velocity(const vector<Cartosphere::Point>& points,
	vector<FL3>& velocities) const
{
	// TODO: Implement
}
