
#include <iostream>
#include <sstream>

#include "cartosphere/mesh.hpp"

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
	unsigned const total = 5;
	for (unsigned i = 0; i <= total; ++i)
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
		unsigned euler = stats.nPoint + stats.nTriangle - stats.nEdge;
		std::cout << "V - E + F = " << stats.nPoint << " - " << stats.nEdge << " + "
			<< stats.nTriangle << " = " << euler << "\n";
		std::cout << "Area ratio = " << stats.areaElementDisparity
			<< " (max " << stats.areaElementMax
			<< ", min " << stats.areaElementMin << ")\n";
	}

	std::cout << "\nLimit = 4*pi = " << (4 * M_PI) << "\n";

	return 0;
}
