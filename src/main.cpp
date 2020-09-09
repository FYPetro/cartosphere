
#include <iostream>

#include "cartosphere/mesh.hpp"

/* **************************** MAIN ENTRY POINT **************************** */
int
main(int argc, char **argv)
{
	if (argc < 2)
	{
		return 0;
	}

	std::cout << "Input file: " << argv[1] << "\n";
	// Investigate Triangle and MFEM
	Cartosphere::TriangularMesh mesh(argv[1]);
	if (!mesh.isLoaded())
	{
		auto messages = mesh.getMessages();
		for (auto msg : messages)
		{
			std::cout << msg << "\n";
		}
		return 0;
	}

	mesh.reportAreas();
	auto messages = mesh.getMessages();
	for (auto msg : messages)
	{
		std::cout << msg << "\n";
	}

	return 0;
}
