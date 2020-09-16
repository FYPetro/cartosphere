
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

	// Iterative refinement
	unsigned total = 6;
	for (unsigned i = 0; i <= total; ++i)
	{
		std::cout << "Refinement Level " << i << " "
			"Euclidean Area = " << mesh.areaEuclidean() << "\n";

		if (i < total)
		{
			mesh.refine();
			std::string name;
			{
				std::stringstream sst;
				sst << argv[1] << ".r" << i;
				name = sst.str();
			}
			mesh.save(name);
		}
		else
		{
			std::cout << "Limit = 4*pi = " << (4 * M_PI) << "\n";
		}
	}

	return 0;
}
