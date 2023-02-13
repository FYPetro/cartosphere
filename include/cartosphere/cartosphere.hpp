
#ifndef __CARTOSPHERE_HPP__
#define __CARTOSPHERE_HPP__

#include <fftw3.h>
typedef double fftw_real;

#include "cartosphere/solver.hpp"

namespace Cartosphere
{
	// Spherical cartogram scheme
	template <typename DerivedType>
	class SolverWrapper
	{
	public:
		// Transform a vector of points, timestepping from the solver
		void transform(vector<Cartosphere::Point>& points)
		{
			int iteration = 0;
			double timeElapsed = 0;

			// Prepare to loop
			double timestep = firstTimestep;
			double epsilon = DoubleMaximum;
			double maxDistance = DoubleMaximum;
			vector<Cartosphere::Point> previousPoints;

			// Loop while conditions unchange
			bool notExpired = true;
			bool notConvergent = true;
			vector<FL3> velocities(points.size());
			while (notExpired && notConvergent)
			{
				// Compute velocity field
				previousPoints = points;
				velocity(points, velocities);

				// Use velocity field to perform time step.
				FL3 travel;
				double travelDistance;
				maxDistance = 0;
				for (size_t i = 0; i < points.size(); ++i)
				{
					travel = timestep * velocities[i];
					points[i].move(travel);

					travelDistance = travel.norm2();
					if (travelDistance > maxDistance)
					{
						maxDistance = travelDistance;
					}
				}

				// Prepare for next iteration
				++iteration;
				advance_solver(timeElapsed, timestep);
				timeElapsed += timestep;
				timestep *= ratioTimestep;

				// Judge loop criterions
				notExpired = true;
				notConvergent = maxDistance > epsDistance;

				// std::cout << "Iteration #" << iteration << "\n"
				// 	<< "\t" << "Time Elapsed " << timeElapsed << "\n"
				// 	<< "\t" << "Max Distance " << maxDistance << "\n";
			}

			// std::cout << "Iteration stopped.\n";
		}

	protected:
		// Initial condition, defaults the zero function
		Cartosphere::Function initFunction =
			[](const Cartosphere::Point& x) { return 0; };

		// Initial timestep size
		double firstTimestep = 1e-2;

		// Timestep size common ratio, defaults to uniform marching
		double ratioTimestep = 1;

		// Convergence criterion: maximum distance
		double epsDistance = 1e-6;

	public:
		// Advance solver
		void advance_solver(double time, double delta)
		{
			(reinterpret_cast<DerivedType*>(this))->advance_solver(time, delta);
		}

		// Advance solver
		void velocity(const vector<Cartosphere::Point>& points, vector<FL3>& velocities) const
		{
			return (reinterpret_cast<const DerivedType*>(this))->velocity(points, velocities);
		}

		// Get/Set func_ic
		Cartosphere::Function get_initial_condition() const { return initFunction; }
		void set_initial_condition(const Cartosphere::Function& f) { initFunction = f; }

		// Get/Set firstTimeStep
		double get_first_timestep() const { return firstTimestep; }
		void set_first_timestep(double t) { if (t > 0) firstTimestep = t; }

		// Get/Set ratioTimeStep
		double get_ratio_timestep() const { return ratioTimestep; }
		void set_ratio_timestep(double r) { if (r > 0) ratioTimestep = r; }

		// Get/Set epsDistance
		double get_eps_distance() const { return epsDistance; }
		void set_eps_distance(double e) { if (e > 0) epsDistance = e; }
	};

	// Spectral cartogram generator
	class SpectralGlobe : public SolverWrapper<Cartosphere::SpectralGlobe>
	{
	public:
		// Default constructor
		SpectralGlobe() : B(0), N(0) {}

	public:
		// Initialize solver
		void initialize_solver();

		// Advance solver
		void advance_solver(double time, double delta);

		// Compute velocity
		void velocity(const vector<Cartosphere::Point>& points,
			vector<FL3>& velocities) const;

	protected:
		// Data at time 0 and time t
		vector<double> init_data;
		vector<double> time_data;
		
		// Fourier at time 0 and time t
		vector<double> init_hats;
		vector<double> time_hats;

		// Temporary allocations for S2 transformations
		double* ws2;
		fftw_real* ipad;
		fftw_plan idct, idst;

		// Gradient field at cell centers (theta_j^*,phi_k^*)
		vector<double> time_dp;
		vector<double> time_da;
		vector<FL3> time_grad;

		// Bandlimit and data size
		int B;
		int N;

	public:
		// Get/Set bandlimit: must be a whole power of 2 and even
		int get_bandlimit() const { return B; }
		void set_bandlimit(int B) { if (B > 0) this->B = B; }
	};

	// Finite element cartogram generator
	class FiniteElementGlobe : public SolverWrapper<Cartosphere::FiniteElementGlobe>
	{
	public:
		// Initialize solver
		void initialize_solver();
		
		// Advance solver
		void advance_solver(double time, double delta);

		// Compute velocity
		void velocity(const vector<Cartosphere::Point>& points,
			vector<FL3>& velocities) const;
	};
}

#endif // !__CARTOSPHERE_HPP__
