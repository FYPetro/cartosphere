
#include <iostream>

#include "cartosphere/mesh.hpp"

/* **************************** MAIN ENTRY POINT **************************** */
int
main(int argc, char **argv)
{
	Cartosphere::Preimage A, B, C;

	A.p = M_PI / 4;
	A.a = 0;

	B.p = M_PI / 4;
	B.a = 2 * M_PI / 3;

	C.p = M_PI / 4;
	C.a = 4 * M_PI / 3;

	Cartosphere::Triangle ABC;

	ABC.A = A;
	ABC.B = B;
	ABC.C = C;

	std::cout << "A = " << ABC.A.image().x << ", " << ABC.A.image().y << ", " << ABC.A.image().z << "\n";
	std::cout << "B = " << ABC.B.image().x << ", " << ABC.B.image().y << ", " << ABC.B.image().z << "\n";
	std::cout << "C = " << ABC.C.image().x << ", " << ABC.C.image().y << ", " << ABC.C.image().z << "\n";

	std::cout << "Area Spherical = " << ABC.area() << "\n";
	std::cout << "Area Euclidean = " << ABC.areaEuclidean() << "\n";

	if (argc > 1)
	{
		std::cout << "Input file:" << argv[1] << "\n";
		// Investigate Triangle and MFEM
		Cartosphere::TriangularMesh mesh(argv[1]);
	}
	
	return 0;
}















