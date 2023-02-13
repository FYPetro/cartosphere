
#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#include "cartosphere/mesh.hpp"

#include <fftw3.h>

namespace Cartosphere
{
	class SteadyStateSolver
	{
	public:
		// Set mesh
		void set(const TriangularMesh& m)
		{
			_mesh = m;
			_solution.assign(_mesh.statistics().V, 0);
		}
		// Get the solution
		vector<double> get() const
		{
			return _solution;
		}
		// Output the linear system
		void debug(const string& name)
		{
			_debug = name;
		}
		// Solve the steady state system
		void solve(Function f);

	protected:
		// Finite-element Mesh
		TriangularMesh _mesh;
		// Solution
		vector<double> _solution;
		// File name debug
		string _debug;
	};

	class TimeDependentSolver
	{
	public:
		// Set mesh
		void set(const TriangularMesh& m)
		{
			_m = m;

			// Build relevant matrices
			_m.fill(_A, _M, Cartosphere::Triangle::Integrator::Refinement5);

			// Attempt to correct the matrix A
			for (int k = 0; k < _A.outerSize(); ++k)
			{
				SparseMatrixRowMajor::InnerIterator it_diag;
				double sum_offdiag = 0;
				for (SparseMatrixRowMajor::InnerIterator it(_A, k); it; ++it)
				{
					it.row();   // row index
					it.col();   // col index (here it is equal to k)

					// Locate the diagonal element or else accumulate
					if (it.row() == it.col())
					{
						it_diag = it;
					}
					else
					{
						sum_offdiag += it.value();
					}
				}
				it_diag.valueRef() = -sum_offdiag;
			}

			_a = ColVector(_A.cols());
			_v = _m.vertices();
		}

		// Set mesh
		void set(Function f)
		{
			_m.fill(_b, f, Cartosphere::Triangle::Integrator::Refinement5);
		}

		// Set mesh
		void initialize(Function f)
		{
			std::transform(_v.begin(), _v.end(), _a.begin(), f);
		}

		// Advance using forward Euler.
		double advance(double timestep)
		{
			SparseMatrixRowMajor LHS = _A + _M / timestep;
			ColVector RHS = _b + _M / timestep * _a;

			SolverBiCGSTAB s(LHS);
			ColVector a = s.solve(RHS);

			ColVector _aprev = _a;
			_a = a;

			// Update the nodal values but convert Vector to vector<double>
			vector<double> a_vec(_a.size(), 0);
			for (int i = 0; i < _a.size(); ++i)
			{
				a_vec[i] = _a[i];
			}
			_m.set(a_vec);

			double dist = 0;
			for (int i = 0; i < _a.size(); ++i)
			{
				double dist_abs = abs(_aprev[i] - _a[i]);
				if (dist_abs > dist)
				{
					dist = dist_abs;
				}
			}
			return dist;
		}

		// Advance using Crank-Nicolson.
		double advanceCrankNicolson(double timestep);

		// Velocity
		vector<FL3> velocity(const vector<Point>& p) const;

	protected:
		// Finite-element Mesh
		TriangularMesh _m;

		// List of vertices
		vector<Point> _v;

		// Matrix for internal calculation
		SparseMatrixRowMajor _A, _M;

		// Vectors for internal calculation
		ColVector _b, _a;
	};
}

#endif // !__SOLVER_HPP__
