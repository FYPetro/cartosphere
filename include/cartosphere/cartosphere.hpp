
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
		// A snapshot of the solver at various timesteps
		struct Snapshot
		{
			// [time_begin, time_final]
			double time_begin;
			double time_final;
			// Length of interval
			double duration;
			// Extreme velocity and distance
			double max_speed;
			double max_distance;
			// Points and velocity fields
			vector<Cartosphere::Point> points;
		};

	public:
		// Transform a vector of points, timestepping from the solver
		void transform(vector<Cartosphere::Point>& points)
		{
			history.clear();
			
			// Start with only the initial position recorded
			Snapshot status;
			{
				status.time_begin = 0;
				status.time_final = 0;
				status.duration = 0;
				status.max_speed = 0;
				status.max_distance = 0;
				if (recordTrajectory)
				{
					status.points = points;
				}
			}
			history.push_back(status);

			// Prepare to loop
			double timeElapsed = 0;
			double timestep = firstTimestep;
			double epsilon = DoubleMaximum;
			double maxDistance = DoubleMaximum;

			// Loop while conditions unchange
			bool isExpired, isConvergent;
			vector<FL3> velocities(points.size());
			for (int iteration = 0; iteration < maxIterations; ++iteration)
			{
				// Compute velocity field
				advance_solver(timeElapsed, 0);
				velocity(points, velocities);

				// Use velocity field to perform time step.
				FL3 travel;
				double travelDistance;
				double slack;
				maxDistance = 0;
				for (size_t i = 0; i < points.size(); ++i)
				{
					slack = (1 - exp(-2 * timestep)) / (2 * timestep);
					travel = timestep * velocities[i] * slack;
					points[i].move(travel);

					travelDistance = travel.norm2();
					if (travelDistance > maxDistance)
					{
						maxDistance = travelDistance;
					}
				}

				// DO OUR OWN THING
				status.time_begin = timeElapsed;
				status.time_final = timeElapsed + timestep;
				status.duration = timestep;
				status.max_speed = maxDistance / timestep;
				status.max_distance = maxDistance;
				if (recordTrajectory)
				{
					status.points = points;
				}
				history.push_back(status);

				// Prepare for next iteration
				timeElapsed += timestep;
				if (timeAdaptivity)
				{
					timestep = firstTimestep * exp(2 * timeElapsed);
				}
				else
				{
					timestep *= ratioTimestep;
				}

				// Judge loop criterions
				{
					isExpired = timeElapsed > 100;
					isConvergent = maxDistance < epsDistance;
				}
				if (isExpired || isConvergent)
				{
					break;
				}
			}
		}

	protected:
		// A list of all snapshots
		vector<Snapshot> history;

		// Initial condition, defaults the zero function
		Cartosphere::Function initFunction = [](const Cartosphere::Point& x) {
			return 0;
		};

		// Record trajectory?
		bool recordTrajectory = false;

		// Adaptively compute time?
		bool timeAdaptivity = false;

		// Initial timestep size
		double firstTimestep = 1e-3;

		// Timestep size common ratio, defaults to uniform marching
		double ratioTimestep = 1;

		// Convergence criterion: maximum distance
		double epsDistance = 1e-6;

		// Convergence criterion: maximum number of iterations
		double maxIterations = std::numeric_limits<int>::max();

	public:
		// Advance solver
		void advance_solver(double time, double delta)
		{
			(reinterpret_cast<DerivedType*>(this))->advance_solver(time, delta);
		}

		// Compute velocity at given points
		void velocity(const vector<Cartosphere::Point>& points, vector<FL3>& velocities) const
		{
			return (reinterpret_cast<const DerivedType*>(this))->velocity(points, velocities);
		}

		// Output a report
		void format_matlab(const string& prefix) const
		{
			ofstream ofs;
			
			// Prepare file names
			string file_name = prefix + ".m";
			string data_name = prefix + "_data.m";

			// Prepare for logging
			Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

			// Prepare the prefix_data.m file
			ofs.open(data_name);
			if (recordTrajectory)
			{
				for (int i = 0; i < history.size(); ++i)
				{
					auto& snapshot = history[i];

					size_t N = snapshot.points.size();
					ofs << "points(:,:," << (i + 1) << ") = [\n";
					for (size_t j = 0; j < N; ++j)
					{
						ofs << " " << snapshot.points[j].x()
							<< " " << snapshot.points[j].y()
							<< " " << snapshot.points[j].z();
						if (j + 1 == N)
						{
							ofs << "]";
						}
						ofs << "; % snapshot " << i << " point " << j << "\n";
					}
				}
			}
			ofs.close();

			// Prepare the prefix.m file
			ofs.open(file_name);
			{
				ofs << "%% Loading data\n"
					<< "clear points\n"
					<< prefix << "_data" << "\n\n";

				ofs << "%% Drawing trajectories\n"
					<< "sphere;axis equal tight\n"
					<< "xlabel('x');ylabel('y');zlabel('z')\n\n"
					<< "hold on\n"
					<< "for i = 1:size(points,1)\n"
					<< "\t" << "trajectory = squeeze(points(i,:,:))';\n"
					<< "\t" << "XYZ = num2cell(trajectory,1);\n"
					<< "\t" << "plot3(XYZ{:});\n"
					<< "end\n"
					<< "hold off\n\n";

				ofs << "%% Summary\n"
					<< "% Iterations\n";
				for (size_t i = 0; i + 1 < history.size(); ++i)
				{
					auto& snapshot = history[i];

					ofs << "t(" << (i + 1) << ", 1:2) = ["
						<< snapshot.time_begin << ", " << snapshot.time_final << "]; \n"
						<< "v(" << (i + 1) << ",1) = "
						<< snapshot.max_speed << ";\n"
						<< "d(" << (i + 1) << ",1) = "
						<< snapshot.max_distance << ";\n";
				}
			}
			ofs.close();
		}

	public:
		// Get/Set func_ic
		Cartosphere::Function get_initial_condition() const { return initFunction; }
		void set_initial_condition(const Cartosphere::Function& f) { initFunction = f; }

		// Enable/Disable snapshot
		void enable_snapshot() { recordTrajectory = true; }
		void disable_snapshot() { recordTrajectory = false; }

		// Enable/Disable adaptivity
		void enable_time_adaptivity() { timeAdaptivity = true; }
		void disable_time_adaptivity() { timeAdaptivity = false; }

		// Get/Set firstTimeStep
		double get_first_timestep() const { return firstTimestep; }
		void set_first_timestep(double t) { if (t > 0) firstTimestep = t; }

		// Get/Set ratioTimeStep
		double get_ratio_timestep() const { return ratioTimestep; }
		void set_ratio_timestep(double r) { if (r > 0) ratioTimestep = r; }

		// Get/Set epsDistance
		double get_eps_distance() const { return epsDistance; }
		void set_eps_distance(double e) { if (e >= 0) epsDistance = e; }

		// Get/Set maxIterations
		double get_max_iterations() const { return maxIterations; }
		void set_max_iterations(int i) { if (i >= 0) maxIterations = i; }
	};

	// Spectral cartogram generator
	class SpectralGlobe : public SolverWrapper<Cartosphere::SpectralGlobe>
	{
	public:
		// Default constructor
		SpectralGlobe() {}

		~SpectralGlobe() { cleanup(); }

	public:
		// Initialize solver
		void initialize_solver();

		// Cleanup dynamic allocated data
		void cleanup();

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

		// Partials and gradient at time t
		vector<double> time_dp;
		vector<double> time_da;

		// Pole data
		double time_data_north = 0;
		double time_data_south = 0;
		FL3 time_grad_north = { 0, 0, 0 };
		FL3 time_grad_south = { 0, 0, 0 };

		// Allocations for S2 transformations
		vector<double> ws2;
		fftw_real* ipad = nullptr;
		fftw_plan idct = NULL;
		fftw_plan idst = NULL;

		// Bandlimit and data size
		int B = 0;
		int N = 0;

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
