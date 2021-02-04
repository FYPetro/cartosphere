
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include "cartosphere/mesh.hpp"

typedef Eigen::SparseMatrix<FLP, Eigen::RowMajor> CSR;
typedef Eigen::Triplet<FLP> Entry;
typedef std::vector<Entry> Entries;
typedef Eigen::Matrix<FLP, Eigen::Dynamic, 1> VEC;
typedef Eigen::BiCGSTAB<CSR> Solver;

void build_system(const Cartosphere::TriangularMesh& mesh, CSR& A, VEC& b)
{
	// Temporarily use this as a testing ground for solving matrix systems
	int N = 100;
	A.resize(N, N);
	b.resize(N);

	Entries e;

	for (int i = 0; i < N; ++i)
		e.emplace_back(i, i, -2);

	for (int i = 0; i < N - 1; ++i)
		e.emplace_back(i, i + 1, 1);

	for (int i = 0; i < N - 1; ++i)
		e.emplace_back(i + 1, i, 1);

	A.setFromTriplets(e.begin(), e.end());
	b[0] = -1;
	for (size_t i = 1; i < N - 1; ++i)
		b[i] = 0;
	b[N - 1] = -1;
}


/* Demo with no arguments */
int demo()
{
	std::string file = "icosahedron.5.csm";

	// Load mesh from file
	Cartosphere::TriangularMesh mesh(file);
	if (mesh.isReady())
	{
		std::cout << "Loaded mesh from file: " << file << "\n\n";

		// Print statistics about the mesh
		auto stats = mesh.statistics();
		size_t euler = stats.nPoint + stats.nTriangle - stats.nEdge;

		std::cout << "Statistics:\n"
			<< "    Euler: V - E + F = " << stats.nPoint << " - " << stats.nEdge
			<< " + " << stats.nTriangle << " = " << euler << "\n"
			<< "    Area ratio: " << stats.areaElementDisparity
			<< " (max " << stats.areaElementMax
			<< ", min " << stats.areaElementMin << ")\n";
	}
	else
	{
		// Print error messages
		for (auto& msg : mesh.getMessages())
		{
			std::cout << msg << "\n";
		}
		return 0;
	}

	CSR A;
	VEC b;
	build_system(mesh, A, b);

	Solver s(A);
	auto x = s.solve(b);

	std::cout << "# Iterations:    " << s.iterations() << std::endl;
	std::cout << "Estimated Error: " << s.error() << std::endl;
	for (int i = 0; i < x.size(); ++i)
		std::cout << x[i] << "\n";

	// Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(A);
	// Eigen::VectorXd x = chol.solve(b);

	return 0;
}

/* **************************** MAIN ENTRY POINT **************************** */
int main(int argc, char** argv)
{
	// The parsing of arguments
	if (argc < 2)
	{
		return demo();
	}

	return 0;
}
