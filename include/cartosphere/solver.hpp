
#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#include "cartosphere/mesh.hpp"
#include "fftw-3.3.5/api/fftw3.h"

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
		// Advance
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
		FLP advanceCrankNicolson(FLP timestep);
		// Velocity
		std::vector<FL3> velocity(const std::vector<Point>& p) const;
		
	protected:
		// Finite-element Mesh
		TriangularMesh _m;
		// List of vertices
		std::vector<Point> _v;
		// Matrix for internal calculation
		Matrix _A, _M;
		// Vectors for internal calculation
		Vector _b, _a;
	};

	// The old spectral solver using S2kit
	class SpectralSolver
	{
	public:
		// Constructor using just a bandlimit
		SpectralSolver() = default;
		// Destructor
		~SpectralSolver();

	public:
		
		void parse(const std::string& path);

		void execute();

		std::string inputSummary() const;

		std::string outputSummary() const;

	protected:
		// Bandlimit
		int _bandlimit;
		// Initialize workspace

		// Destroy workspace

	private:
		// Bandlimit allocated
		int _B = 0;
		// Memory for diffusion solver
		std::unique_ptr<FLP> _u0, _ut;
		// Memory for spherical harmonic transformations
		std::unique_ptr<FLP> _re0, _im0, _ret, _imt;
		// Memory reserved for S2kit
		std::unique_ptr<FLP> _ws, _wt;
		// FFTW Plan for Fast Fourier Transform
		fftw_plan _fft = nullptr;
		// FFTW Plan for Inverse Fast Fourier Transform
		fftw_plan _ifft = nullptr;
		// FFTW Plan for DCT-II
		fftw_plan _dct = nullptr;
		// FFTW Plan for DCT-III
		fftw_plan _idct = nullptr;
	};
}

#endif // !__SOLVER_HPP__
