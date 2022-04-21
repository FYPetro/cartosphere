
#include "cartosphere/solver.hpp"

#include <fstream>

using Cartosphere::TimeDependentSolver;
using Cartosphere::Arc;
using Cartosphere::Point;

void
Cartosphere::SteadyStateSolver::solve(Function f)
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
