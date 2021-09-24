
#include "cartosphere/solver.hpp"

#include <fstream>

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