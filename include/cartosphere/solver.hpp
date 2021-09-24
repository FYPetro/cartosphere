
#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#include "cartosphere/mesh.hpp"

namespace Cartosphere
{
	class SteadyStateSolver
	{
	public:
		// Set mesh
		void set(const TriangularMesh& m)
		{
			_mesh = m;
			_solution.assign(_mesh.statistics().V, FLP(0));
		}
		// Get the solution
		std::vector<FLP> get() const
		{
			return _solution;
		}
		// Output the linear system
		void debug(const std::string& name)
		{
			_debug = name;
		}
		// Solve the steady state system
		void solve(Function f);

	protected:
		// Finite-element Mesh
		TriangularMesh _mesh;
		// Solution
		std::vector<FLP> _solution;
		// File name debug
		std::string _debug;
	};
}

#endif // !__SOLVER_HPP__
