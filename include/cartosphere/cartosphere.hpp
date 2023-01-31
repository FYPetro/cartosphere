
#ifndef __CARTOSPHERE_HPP__
#define __CARTOSPHERE_HPP__

#include "cartosphere/solver.hpp"

namespace Cartosphere
{
	// Spherical cartogram scheme
	template <typename _Solver>
	class SolverWrapper
	{
	public:
		// Transform a vector of points, timestepping from the solver
		void transform(vector<Cartosphere::Point>& points)
		{
			int iteration = 0;
			FLP timeElapsed = 0;

			// Prepare to loop
			FLP timestep = firstTimestep;
			FLP epsilon = std::numeric_limits<FLP>::max();
			FLP maxDistance = std::numeric_limits<FLP>::max();
			vector<Cartosphere::Point> previousPoints;

			// Loop while conditions unchange
			bool notExpired = true;
			bool notDivergent = true;
			while (notExpired && notDivergent)
			{
				// Compute velocity field
				previousPoints = points;
				vector<FL3> velocity = solver.velocity(points);

				// Use velocity field to perform time step.
				FL3 travel;
				FLP travelDistance;
				maxDistance = 0;
				for (size_t i = 0; i < points.size(); ++i)
				{
					travel = timestep * velocity[i];
					points[i].move(travel);

					travelDistance = travel.norm2();
					if (travelDistance > maxDistance)
					{
						maxDistance = travelDistance;
					}
				}

				// Prepare for next iteration
				++iteration;
				solver.advance(timestep);
				timeElapsed += timestep;
				timestep *= ratioTimestep;

				// Judge loop criterions
				notExpired = true;
				notDivergent = maxDistance > epsDistance;

				std::cout << "Iteration #" << iteration << "\n"
					<< "\t" << "Time Elapsed " << timeElapsed << "\n"
					<< "\t" << "Max Distance " << maxDistance << "\n";
			}

			std::cout << "Iteration stopped.\n";
		}

	protected:
		// Solver instance
		_Solver solver;

		// Initial timestep size
		FLP firstTimestep;

		// Timestep common ratio
		FLP ratioTimestep;

		// Convergence criterion: maximum distance
		FLP epsDistance;

	public:
		// Initialize solver
		void initialize_solver() const;

		// Get/Set firstTimeStep
		FLP get_first_timestep() const { return firstTimestep; }
		void set_first_timestep(FLP timestep) { firstTimestep = timestep; }

		// Get/Set ratioTimeSTep
		FLP get_ratio_timestep() const { return ratioTimestep; }
		void set_ratio_timestep(FLP ratio) { if (ratio > 0) ratioTimestep = ratio; }

		// Get/Set epsDistance
		FLP get_eps_distance() const { return epsDistance; }
		void set_eps_distance(FLP eps) { if (eps > 0) epsDistance = eps; }
	};

	// Spectral cartogram generator
	class SpectralGlobe : public SolverWrapper<Cartosphere::SpectralSolver>
	{

	};

	// Finite element cartogram generator
	class FiniteElementGlobe : public SolverWrapper<Cartosphere::TimeDependentSolver>
	{

	};
}

#endif // !__CARTOSPHERE_HPP__
