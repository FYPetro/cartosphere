
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
				Matrix::InnerIterator it_diag;
				FLP sum_offdiag = 0;
				for (Matrix::InnerIterator it(_A, k); it; ++it)
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

			_a = Vector(_A.cols());
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
		FLP advance(FLP timestep)
		{
			Matrix LHS = _A + _M / timestep;
			Vector RHS = _b + _M / timestep * _a;

			Solver s(LHS);
			Vector a = s.solve(RHS);

			Vector _aprev = _a;
			_a = a;

			// Update the nodal values but convert Vector to std::vector<FLP>
			std::vector<FLP> a_vec(_a.size(), 0);
			for (int i = 0; i < _a.size(); ++i)
			{
				a_vec[i] = _a[i];
			}
			_m.set(a_vec);

			FLP dist = 0;
			for (int i = 0; i < _a.size(); ++i)
			{
				FLP dist_abs = abs(_aprev[i] - _a[i]);
				if (dist_abs > dist)
				{
					dist = dist_abs;
				}
			}
			return dist;
		}

		// Advance using Crank-Nicolson.
		FLP advanceCrankNicolson(FLP timestep);

		// Velocity
		vector<FL3> velocity(const std::vector<Point>& p) const;

	protected:
		// Finite-element Mesh
		TriangularMesh _m;

		// List of vertices
		vector<Point> _v;

		// Matrix for internal calculation
		Matrix _A, _M;

		// Vectors for internal calculation
		Vector _b, _a;
	};

	// The old spectral solver using S2kit
	class SpectralSolver
	{
	public:
		// Constructor using default bandlimit
		SpectralSolver() : SpectralSolver(0) {};

		// Constructor using custom bandlimit
		SpectralSolver(int B) : bandlimit(B) { ws_initialize(); }

		// Disable copy constructor
		SpectralSolver(const SpectralSolver &other) : bandlimit(solver.bandlimit)
		{
			int offset = 4 * pow(bandlimit, 2);
			std::copy(other.u0, other.u0 + offset, u0);
			std::copy(other.ut, other.ut + offset, ut);
			std::copy(other.r0, other.r0 + offset, r0);
			std::copy(other.rt, other.rt + offset, rt);

			offset = 10 * pow(bandlimit, 2) + 24 * bandlimit;
			std::copy(other.ws, other.ws + offset, ws);

			offset = 4 * bandlimit;
			std::copy(other.wt, other.wt + offset, wt);
		}

		// Destructor
		~SpectralSolver() { if (ws != nullptr) ws_destroy(); }

	public:
		// Velocity field calculation
		vector<FL3> velocity(vector<Cartosphere::Point> &points) const;

	protected:
		// Bandlimit
		int bandlimit;

		// Memory for diffusion solver
		std::unique_ptr<FLP> u0, ut;

		// Memory for Fourier transformations
		std::unique_ptr<FLP> r0, rt;

		// Memory reserved for S2kit
		std::unique_ptr<FLP> ws, wt;

		// FFTW plans
		fftw_plan dct, idct, ws_dct, ws_idct;

		// Initialize workspace
		void ws_initialize();

		// Execute workspace
		void ws_execute();

		// Destroy workspace
		void ws_destroy();

	public:
		// Get/Set bandlimit
		int get_bandlimit() const { return bandlimit; }
		void set_bandlimit(int B) { if (B > 0) bandlimit = B; }
	};
}

#endif // !__SOLVER_HPP__
