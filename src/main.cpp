
#include <iostream>
#include <sstream>

#include "cartosphere/mesh.hpp"

// Tests the matrix
#include "cartosphere/algebra.hpp"

void test_matrix()
{
	Cartosphere::Matrix A(5, 3);
	A = { 1, 2, 3, 4, 5, 1, 3, 5, 2, 4, 1, 4, 2, 5, 3 };

	std::vector<FLP> b = { -10, 12, 14, 16, 18 };
	std::vector<FLP> x;
	A.solve(x, b);

	for (size_t i = 0; i < x.size(); ++i)
	{
		std::cout << x[i] << "\n";
	}
}

/* **************************** MAIN ENTRY POINT **************************** */
int
main(int argc, char **argv)
{
	// The parsing of arguments
	if (argc < 2)
	{
		return 0;
	}

	std::cout << "Input file: " << argv[1] << "\n";
	Cartosphere::TriangularMesh mesh(argv[1]);
	if (!mesh.isReady())
	{
		for (auto &msg : mesh.getMessages())
		{
			std::cout << msg << "\n";
		}
		return 0;
	}

	mesh.format(std::string(argv[1]) + ".obj");

	// Iterative refinement
	size_t const total = 3;
	for (size_t i = 0; i <= total; ++i)
	{
		std::cout << "\nRefinement Level " << i << "\n"
			"Euclidean Area = " << mesh.areaEuclidean() << "\n"
			"Spherical Area = " << mesh.area() << "\n";

		if (i > 0)
		{
			mesh.refine();

			std::string name;
			{
				std::stringstream sst;
				sst << argv[1] << ".r" << i;
				name = sst.str();
			}
			mesh.format(name + ".obj");

			mesh.save(name);
		}

		// Print statistics about the mesh
		auto stats = mesh.statistics();
		size_t euler = stats.nPoint + stats.nTriangle - stats.nEdge;
		std::cout << "V - E + F = " << stats.nPoint << " - " << stats.nEdge
			<< " + " << stats.nTriangle << " = " << euler << "\n";
		std::cout << "Area ratio = " << stats.areaElementDisparity
			<< " (max " << stats.areaElementMax
			<< ", min " << stats.areaElementMin << ")\n";
	}

	std::cout << "\nLimit = 4*pi = " << (4 * M_PI) << "\n";

	test_matrix();

	return 0;
}
