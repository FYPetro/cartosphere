
#include "cartosphere/cartosphere.hpp"
using Cartosphere::Point;
using Cartosphere::SpectralGlobe;
using Cartosphere::FiniteElementGlobe;

#include "cartosphere/dsht.hpp"

void
SpectralGlobe::initialize_solver()
{
	// If the non-zero bandlimit is new
	if (N != 2 * B && B > 0)
	{
		N = 2 * B;
		// Resize existing
		init_data.resize(N * N);
		time_data.resize(N * N);
		init_hats.resize(B * B);
		time_hats.resize(B * B);
		time_dp.resize(N * N);
		time_da.resize(N * N);
		time_grad.resize(N * N);
		// Free old dynamic allocations
		delete[] ws2;
		fftw_free(ipad);
		fftw_destroy_plan(idct); fftw_destroy_plan(idst);
		// Allocate new stuff
		ws2 = cs_make_ws2(B);
		ipad = fftw_alloc_real(N * N * 2);
		cs_ids2ht_plans(B, ipad, &idct, &idst);
	}

	// Clear data and hats
	if (N > 0)
	{
		// Data and hats storage
		memset(init_data.data(), 0, N * N * sizeof(double));
		memset(time_data.data(), 0, N * N * sizeof(double));
		memset(init_hats.data(), 0, B * B * sizeof(double));
		memset(time_hats.data(), 0, B * B * sizeof(double));
		
		// Sample initial condition
		double phi, theta;
		for (int j = 0; j < N; ++j)
		{
			theta = M_PI / N * (j + .5);
			for (int k = 0; k < N; ++k)
			{
				phi = M_PI / B * (k + .5);
				Point P(theta, phi);
				init_data[N * j + k] = initFunction(P);
			}
		}

		// Compute initial Fourier coefficients
		cs_fds2ht(B, init_data.data(), init_hats.data(), ws2);
		// Update initial data and gradient
		advance_solver(0, 0);
	}
}

void
SpectralGlobe::advance_solver(double time, double delta)
{
	// Interval: [0, t]
	double t = time + delta;

	// Compute decayed coefficients
	{
		int l, m, i;
		for (l = 0; l < B; ++l)
		{
			double omega = l * (l + 1);
			for (m = -l; m <= l; ++m)
			{
				i = cs_index2(B, l, m);
				time_hats[i] = init_hats[i] * exp((-omega) * t);
			}
		}
	}

	// Compute homogenized data
	cs_ids2ht(B, time_hats.data(), time_data.data(), ws2, ipad, idct, idst);

	// Compute a velocity field at each grid cell corner
	cs_ids2ht_dp(B, time_hats.data(), time_dp.data(), ws2, ipad, idct, idst);
	cs_ids2ht_da(B, time_hats.data(), time_da.data(), ws2, ipad, idct, idst);
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
				grad.x = time_dp[i] * cos_theta * cos_phi - time_da[i] * sin_phi;
				grad.y = time_dp[i] * cos_theta * sin_phi - time_da[i] * cos_phi;
				grad.z = time_dp[i] * (-sin_theta);
			}
		}
	}
}

void
SpectralGlobe::velocity(const vector<Point>& points, vector<FL3>& velocities) const
{
	// Compute data and velocities at the poles
	double data_north = 0;
	double data_south = 0;
	for (int l = 0; l < B; ++l)
	{
		double q_l = sqrt((2 * l + 1) / (4 * M_PI));
		double q_hat_l = q_l * time_hats[cs_index2(B, l, 0)];
		data_north += q_hat_l;
		data_south += ((l % 2) ? -1 : 1) * q_hat_l;
	}
	FL3 grad_north{ 0, 0, 0 };
	FL3 grad_south{ 0, 0, 0 };
	{
		double phi, diff;
		for (int k = 0; k < N; ++k)
		{
			phi = M_PI / B * (k + 0.5);
			diff = (time_data[k] - data_north);
			grad_north.x += diff * cos(phi);
			grad_north.y += diff * sin(phi);
			diff = (time_data[k] - data_south);
			grad_south.x += diff * cos(phi);
			grad_south.y += diff * sin(phi);
		}
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

		// Obtain data and gradient on the northern side
		double data_n;
		FL3 grad_n;
		if (j_n == -1)
		{
			// The northern edge is degenerate; it is the north pole
			data_n = data_north;
			grad_n = grad_north;
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
			data_s = data_south;
			grad_s = grad_south;
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
		velocities[i] = grad * (-1 / data);
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
