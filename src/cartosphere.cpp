
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
	for (int l = 0; l < B; ++l)
	{
		int eigenvalue = -l * (l + 1);
		for (int m = -l; m <= l; ++m)
		{
			int i = cs_index2(B, l, m);
			H[i] = init_hats[i] * exp(eigenvalue * t);
		}
	}

	// Compute homogenized data
	cs_ids2ht(B, H, D, W, ipad, idct, idst);

	// Compute a velocity field at each grid cell corner
	cs_ids2ht_dp(B, H, P[0], W, ipad, idct, idst);
	cs_ids2ht_da(B, H, P[1], W, ipad, idct, idst);

	// Compute data and velocities at the poles
	time_data_north = 0;
	time_data_south = 0;
	for (int l = 0; l < B; ++l)
	{
		double q_l = sqrt((l + 0.5) / (2 * M_PI));
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
			<< time_grad_north.y << ", " << time_grad_north.z << "]" << "\n"
			<< "  SOUTH: [" << time_grad_south.x << ", "
			<< time_grad_south.y << ", " << time_grad_south.z << "]";
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
		int NW, NE, SW, SE;

		// Compute remainder coordinates for later bilinear interpolation
		j_frac -= j_n;
		k_frac -= k_w;

		double data_north, sin_north, dp_north, da_north;
		// Northern spherical cap
		if (j_n == -1)
		{
			// Mark using magic numbers that the northern edge is degenerate
			NW = -1;
			NE = -1;
			sin_north = 0;
			// Scale j_frac from [0.5,1] to [0,1]
			j_frac = 2 * j_frac - 1;
			// Pick data from the north pole
			data_north = time_data_north;
			// Convert the gradient at the north pole into local coordinates
			double azimuth = M_PI / B * k_e;
			// Compute component along unit tangent at the north pole
			FL3 basis_theta = {cos(azimuth), sin(azimuth), 0};
			dp_north = dot(time_grad_north, basis_theta);
			// Compute component along unit normal at the north pole
			FL3 basis_phi = {-sin(azimuth), cos(azimuth), 0};
			da_north = dot(time_grad_north, basis_phi);
		}
		// Northern edge is non-degenerate
		else
		{
			// Compute the index of the northwest and northeast nodes
			NW = N * j_n + k_w;
			NE = N * j_n + k_e;
			sin_north = sin(M_PI / N * (j_n + 0.5));
			// Perform linear interpolation along the northern edge
			data_north = (1 - k_frac) * time_data[NW] + k_frac * time_data[NE];
			dp_north = (1 - k_frac) * time_dp[NW] + k_frac * time_dp[NE];
			da_north = (1 - k_frac) * time_da[NW] + k_frac * time_da[NE];
		}

		double data_south, sin_south, dp_south, da_south;
		// Southern spherical cap
		if (j_s == N)
		{
			// Mark using magic numbers that the southern edge is degenerate
			SW = -1;
			SE = -1;
			sin_south = 0;
			// Scale j_frac from [0,0.5] to [0,1]
			j_frac = 2 * j_frac;
			// Pick data from the south pole
			data_south = time_data_south;
			// Convert the gradient at the south pole into local coordinates
			double azimuth = M_PI / B * k_e;
			// Compute component along unit tangent at the south pole
			FL3 basis_theta = {-cos(azimuth), -sin(azimuth), 0};
			dp_south = dot(time_grad_south, basis_theta);
			// Compute component along unit normal at the south pole
			FL3 basis_phi = {-sin(azimuth), cos(azimuth), 0};
			da_south = dot(time_grad_south, basis_phi);
		}
		// Southern edge is non-degenerate
		else
		{
			// Compute the index of the southwest and southeast nodes
			SW = N * j_s + k_w;
			SE = N * j_s + k_e;
			sin_south = sin(M_PI / N * (j_s + 0.5));
			// Perform linear interpolation along the southern edge
			data_south = (1 - k_frac) * time_data[SW] + k_frac * time_data[SE];
			dp_south = (1 - k_frac) * time_dp[SW] + k_frac * time_dp[SE];
			da_south = (1 - k_frac) * time_da[SW] + k_frac * time_da[SE];
		}

		// Complete the bilinear interpolation for data and gradient components
		// along the local basis
		double data = (1 - j_frac) * data_north + j_frac * data_south;
		double u = (1 - j_frac) * dp_north + j_frac * dp_south;
		double v = 0;
		if (sin_north == 0)
		{
			v += (1 - j_frac) * da_north;
		}
		else
		{
			v += (1 - j_frac) * da_north / sin_north;
		}
		if (sin_south == 0)
		{
			v += j_frac * da_south;
		}
		else
		{
			v += j_frac * da_south / sin_south;
		}

		// Turn u e_theta + v e_phi into cartesian coordinates
		FL3 grad;
		{
			double cos_theta = cos(P.p());
			double sin_theta = sin(P.p());
			double cos_phi = cos(P.a());
			double sin_phi = sin(P.a());
			// Convert gradient!
			grad.x = u * cos_theta * cos_phi - v * sin_phi;
			grad.y = u * cos_theta * sin_phi + v * cos_phi;
			grad.z = u * (-sin_theta);
		}
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
