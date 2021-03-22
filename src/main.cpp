
#include <iostream>
#include <fstream>

#include "cartosphere/mesh.hpp"

/* Build a linear system for the demo program */
void build_system(const Cartosphere::TriangularMesh& mesh, Matrix& A, Vector& b)
{
	mesh.fill(A);

	auto f = [](const Cartosphere::Point &p)->FLP {
		return p.x() + p.y() + p.z();
	};
	mesh.fill(b, f);
}

/* Demo with no arguments */
int demo()
{
	std::string file = "icosahedron.csm.2.csm";

	// Load mesh from file
	Cartosphere::TriangularMesh mesh(file);
	if (mesh.isReady())
	{
		std::cout << "Loaded mesh from file: " << file << "\n\n";

		// Print statistics about the mesh
		auto stats = mesh.statistics();
		size_t euler = stats.V + stats.F - stats.E;

		std::cout << "Statistics:\n"
			<< "    Euler: V - E + F = " << stats.V << " - " << stats.E
			<< " + " << stats.F << " = " << euler << "\n"
			<< "    Area ratio: " << stats.areaElementDisparity
			<< " (max " << stats.areaElementMax
			<< ", min " << stats.areaElementMin << ")" << std::endl;
	}
	else
	{
		// Print error messages
		for (auto& msg : mesh.getMessages())
		{
			std::cout << msg << "\n";
		}
		std::cout << std::flush;
		return 0;
	}

	// Build and solve Ax=b
	Matrix A;
	Vector b;
	build_system(mesh, A, b);
	Solver s(A);
	auto x = s.solve(b);

	
	if (A.rows() < 0)
	{
		// Print A
		std::cout << "\n";
		{
			Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n   ", "", "", "A = [", "];");
			std::cout << Eigen::MatrixXd(A).format(fmt) << std::endl;
		}
		// Print b
		std::cout << "\n";
		{
			Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n   ", "", "", "b = [", "];");
			std::cout << Eigen::MatrixXd(b).format(fmt) << std::endl;
		}
		// Print x
		std::cout << "\n";
		{
			Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n   ", "", "", "x = [", "];");
			std::cout << Eigen::MatrixXd(x).format(fmt) << std::endl;
		}
	}
	else
	{
		std::ofstream ofs;
		ofs.open("temp.m");
		// Print A
		ofs << "\n";
		{
			Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n   ", "", "", "A = [", "];");
			ofs << Eigen::MatrixXd(A).format(fmt) << std::endl;
		}
		// Print b
		ofs << "\n";
		{
			Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n   ", "", "", "b = [", "];");
			ofs << Eigen::MatrixXd(b).format(fmt) << std::endl;
		}
		// Print x
		ofs << "\n";
		{
			Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n   ", "", "", "x = [", "];");
			ofs << Eigen::MatrixXd(x).format(fmt) << std::endl;
		}
		std::cout << "System too large, wrote to file \"temp.m\"" << std::endl;
	}

	std::cout << "\n";
	std::cout << "Solver Statistics:\n";
	std::cout << "# Iterations:    " << s.iterations() << std::endl;
	std::cout << "Estimated Error: " << s.error() << std::endl;

	return 0;
}

/* Quadrature rule demo */
int demo_quadrature()
{
	// Demo triangle
	Cartosphere::Point P(Cartosphere::Image(0, 0, 1));
	Cartosphere::Point A(Cartosphere::Image(1, 0, 0));
	Cartosphere::Point B(Cartosphere::Image(0, 1, 0));
	Cartosphere::Triangle t(P, A, B);
	Cartosphere::TriangularMesh original(t);

	// Demo function
	Cartosphere::Function f = [](const Cartosphere::Point& p)->FLP {
		return p.x();
	};

	// Print statistics about the mesh
	auto stats = original.statistics();
	size_t euler = stats.V + stats.F - stats.E;
	std::cout << "Statistics:\n"
		<< "    Euler: V - E + F = " << stats.V << " - " << stats.E
		<< " + " << stats.F << " = " << euler << "\n"
		<< "    Area ratio: " << stats.areaElementDisparity
		<< " (max " << stats.areaElementMax
		<< ", min " << stats.areaElementMin << ")\n";

	// Print exact value of integral
	std::cout << "Exact value for integral: " << M_PI_4 << std::endl;

	// Algorithm 1: power-2 refinement of centroid rule
	auto mesh = original;
	std::cout << "Power-2 refinement:\n";
	for (unsigned k = 1; k <= 10; ++k)
	{
		mesh.refine();
		std::cout << "Level " << k << " integral: "
			<< mesh.integrate(f,
				Cartosphere::TriangularMesh::Quadrature::AreaWeighted,
				Cartosphere::Triangle::Integrator::Centroid)
			<< "\n";
	}

	// Algorithm 2: power-2 refinement of three-vertices rule
	mesh = original;
	std::cout << "Power-2 refinement:\n";
	for (unsigned k = 1; k <= 10; ++k)
	{
		mesh.refine();
		std::cout << "Level " << k << " integral: "
			<< mesh.integrate(f,
				Cartosphere::TriangularMesh::Quadrature::AreaWeighted,
				Cartosphere::Triangle::Integrator::ThreeVertices)
			<< "\n";
	}

	return 0;
}

/* Generate weights */
int precompute_weights(std::string path)
{
	Cartosphere::TriangularMesh mesh(path);
	// TODO
	return 0;
}

int refine(std::string path)
{
	Cartosphere::TriangularMesh m(path);
	for (unsigned k = 1; k <= 5; ++k)
	{
		std::stringstream sst;
		sst << path << "." << k << ".csm";
		std::string name = sst.str();
		m.refine();
		m.save(name);
	}
	return true;
}

/* **************************** MAIN ENTRY POINT **************************** */
int main(int argc, char** argv)
{
	// The parsing of arguments
	if (argc < 2)
	{
		return demo();
	}

	// Switch based on the first keyword
	std::string item(argv[1]);
	if (item == "quadrature")
	{
		return demo_quadrature();
	}
	else if (item == "precompute")
	{
		return precompute_weights(argv[2]);
	}
	else if (item == "refine")
	{
		return refine(argv[2]);
	}

	return 0;
}
