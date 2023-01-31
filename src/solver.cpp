
#include "cartosphere/solver.hpp"

#include <fstream>

using namespace Cartosphere;

void
SteadyStateSolver::solve(Function f)
{
	// ***********************
	// Build the linear system
	// ***********************
	
	Matrix A;
	_mesh.fill(A);

	Vector F;
	_mesh.fill(F, f);

	// ***********************
	// Solve the linear system
	// ***********************

	Solver kernel(A);
	Vector x = kernel.solve(F);

	// ******************
	// Process the output
	// ******************
	std::transform(x.begin(), x.end(), _solution.begin(),
		[](auto& v) -> FLP { return v; }
	);

	// DEBUG
	if (!_debug.empty())
	{
		std::ofstream ofs(_debug);
		// Print A
		ofs << "\n";
		{
			Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n   ", "", "", "A = [", "];");
			ofs << Eigen::MatrixXd(A).format(fmt) << std::endl;
		}
		// Print b
		ofs << "\n";
		{
			Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n   ", "", "", "F = [", "];");
			ofs << Eigen::MatrixXd(F).format(fmt) << std::endl;
		}
		// Print x
		ofs << "\n";
		{
			Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n   ", "", "", "x = [", "];");
			ofs << Eigen::MatrixXd(x).format(fmt) << std::endl;
		}
		// Print smallest eigenvalue
		// Print sum of F
		ofs << "\n"
			<< "min(eig(A))\n"
			<< "sum(F)\n";
	}
}

FLP
TimeDependentSolver::advanceCrankNicolson(FLP timestep)
{
	Matrix LHS = _A / 2 + _M / timestep;
	Vector RHS = _b + (_M / timestep - _A / 2) * _a;

	Solver s(LHS);
	Vector a = s.solve(RHS);

	_a = a;

	// Update the nodal values but convert Vector to std::vector<FLP>
	std::vector<FLP> a_vec(_a.size(), 0);
	for (int i = 0; i < _a.size(); ++i)
	{
		a_vec[i] = _a[i];
	}
	_m.set(a_vec);

	return 0;
}

std::vector<FL3>
TimeDependentSolver::velocity(const std::vector<Point>& p) const
{
	std::vector<FL3> u;
	u.reserve(p.size());
	for (size_t k = 0; k < p.size(); ++k)
	{
		FL3 gradient = _m.gradient(p[k]);
		FLP value = _m.interpolate(p[k]);
		u.push_back(-gradient / value);
	}
	return u;

	// for (auto& x : p)
	// {
	//	u.push_back(-_m.gradient(x)/_m.interpolate(x));
	// }
	// return u;
	// std::vector<FL3> u;
	//
	// // Interpolation should happen here
	// for (int i = 0; i < p.size(); ++i)
	// {
	//	FL3 coeff;
	//	UI3 index = _m.indexTriangleVertices(_m.barycentric(p, &coeff));
	//	coeff.normalize();
	//
	//	FL3 ta = Cartosphere::transport(v[index.a], _v[index.a], p);
	//	FL3 tb = Cartosphere::transport(v[index.b], _v[index.b], p);
	//	FL3 tc = Cartosphere::transport(v[index.c], _v[index.c], p);
	//
	//	u[i] = ta * coeff.x + tb * coeff.y + tc * coeff.z;
	// }
	//
	// return u;
}

SpectralSolver::~SpectralSolver()
{

}

void
SpectralSolver::parse(const std::string& path)
{
	// Set a temporary bandlimit to test the code
	_B = 16;
}

void
SpectralSolver::execute()
{
	/*
	// (Re)allocate memory if necessary
	if (_bandlimit != _B)
	{
		// Update stored bandlimit
		_B = _bandlimit;

		// Sampling rate is doubled
		int size = 2 * _B;

		// Initial distribution (2B x 2B)
		_u0.reset(new FLP[size * size]);
		// Diffused distribution (2B x 2B)
		// [used towards fftw plan]
		_ut.reset(fftw_alloc_real(size * size));
		// Real parts of initial harmonic coefficients
		_re0.reset(new FLP[size * size]);
		// Imaginary parts of initial harmonic coefficients
		_im0.reset(new FLP[size * size]);
		// Real parts of diffused harmonic coefficients
		// [used towards fftw plan]
		_ret.reset(fftw_alloc_real(_B * _B));
		// Imaginary parts of diffused harmonic coefficients
		// [used towards fftw plan]
		_imt.reset(fftw_alloc_real(_B * _B));

		// Initialize workspace required by S2kit
		_ws.reset(fftw_alloc_real(10 * _B * _B + 24 * _B));
		// Initialize weights required by S2kit
		_wt.reset(fftw_alloc_real(4 * _B));
		// Fill array with weights
		makeweights(_B, _ws.get());

		// Initialize the fftw plans
		int rank, howmany_rank;
		fftw_iodim dims[1], howmany_dims[1];
		
		// Starting with the FFT plan
		// https://www.fftw.org/fftw3_doc/Guru-Real_002ddata-DFTs.html
		rank = 1;
		dims[0].n = size;
		dims[0].is = 1;
		dims[0].os = size;

		howmany_rank = 1;
		howmany_dims[0].n = size;
		howmany_dims[0].is = size;
		howmany_dims[0].os = 1;

		_fft = fftw_plan_guru_split_dft_r2c(
			rank, dims, howmany_rank, howmany_dims,
			_ut.get(), _ws.get(), _ws.get() + size * size,
			FFTW_MEASURE
		);

		// Then plan the inverse FFT
		// https://www.fftw.org/fftw3_doc/Guru-Real_002ddata-DFTs.html
		rank = 1;
		dims[0].n = size;
		dims[0].is = size;
		dims[0].os = 1;

		howmany_rank = 1;
		howmany_dims[0].n = size;
		howmany_dims[0].is = 1;
		howmany_dims[0].os = size;

		_ifft = fftw_plan_guru_split_dft_c2r(
			rank, dims, howmany_rank, howmany_dims,
			_ws.get(), _ws.get() + size * size, _ut.get(),
			FFTW_MEASURE
		);

		// Plan the DCT
		_dct = fftw_plan_r2r_1d(
			size, _wt.get(), _ut.get(),
			FFTW_REDFT10, FFTW_MEASURE
		);

		// Plan the inverse DCT
		_idct = fftw_plan_r2r_1d(
			size, _wt.get(), _ut.get(),
			FFTW_REDFT01, FFTW_MEASURE
		);
	}
	*/
}

std::string
SpectralSolver::inputSummary() const
{
	std::stringstream sst;

	return sst.str();
}

std::string
SpectralSolver::outputSummary() const
{
	std::stringstream sst;

	return sst.str();
}
