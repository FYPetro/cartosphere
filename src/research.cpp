
#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>

#include "cartosphere/research.hpp"
#include "cartosphere/functions.hpp"
#include "cartosphere/globe.hpp"
#include "cartosphere/shapefile.hpp"

using Cartosphere::Point;
using Cartosphere::TriangularMesh;
using Cartosphere::TimeDependentSolver;
using Cartosphere::Function;

void build_system(const Cartosphere::TriangularMesh& mesh, Matrix& A, Vector& b)
{
	mesh.fill(A);

	auto f = [](const Cartosphere::Point& p)->FLP {
		return p.x() + p.y() + p.z();
	};
	mesh.fill(b, f);
}

int demo()
{
	std::string file = "icosahedron.csm";

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
	Vector x = s.solve(b);


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

	std::vector<FLP> solution;
	solution.reserve(x.size());
	for (int k = 0; k < x.size(); ++k)
	{
		solution.push_back(x[k]);
	}
	mesh.format("demo.obj", solution);

	return 0;
}

int demo_diffusion()
{
	std::string file = "icosahedron.csm";

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

	// Build relevant matrices
	Matrix A, M;
	mesh.fill(A, M, Cartosphere::Triangle::Integrator::Refinement5);

	// Attempt to correct the matrix A
	for (int k = 0; k < A.outerSize(); ++k)
	{
		Matrix::InnerIterator it_diag;
		FLP sum_offdiag = 0;
		for (Matrix::InnerIterator it(A, k); it; ++it)
		{
			it.row();   // row index
			it.col();   // col index (here it is equal to k)

			// Locate the diagonal element or else accumulate
			if (it.row() == it.col())
			{
				it_diag = it;
			}
			else
			{
				sum_offdiag += it.value();
			}
		}
		it_diag.valueRef() = -sum_offdiag;
	}

	auto f = [](const Cartosphere::Point& p)->FLP { return 0; };
	Vector b;
	mesh.fill(b, f, Cartosphere::Triangle::Integrator::Refinement5);

	// Design an initial condition
	Vector v_init(A.cols());
	auto vs = mesh.vertices();

	Cartosphere::Function v_init_f = [](const Cartosphere::Point& p)->FLP {
		return 2 + p.z();
	};
	std::transform(vs.begin(), vs.end(), v_init.begin(), v_init_f);

	// Start!
	FLP tolerance = 1e-6;
	FLP indicator = 1;
	int time_steps = 20;
	FLP time_elapsed = 20;

	int iteration = 0;
	Vector v_prev = v_init, v_curr;

	std::cout << b << "\n";
	std::cout << A << "\n";

	for (int step = 0; step < time_steps; ++step)
	{
		FLP duration = time_elapsed / time_steps;

		Matrix LHS = A + M / duration;
		Vector RHS = b + M / duration * v_prev;

		Solver s(LHS);
		v_curr = s.solve(RHS);

		// Convergence criterion
		indicator = (v_curr - v_prev).norm();

		// Another indicator for debugging only
		FLP max = *std::max_element(std::begin(v_curr), std::end(v_curr)),
			min = *std::min_element(std::begin(v_curr), std::end(v_curr));
		FLP range = max - min;

		std::cout
			<< "Iteration #" << iteration++ << "\n"
			<< " Time Step  " << duration << "\n"
			<< " Indicator " << indicator << "\n"
			<< " Range     " << range << " [" << min << "," << max << "]\n";

		v_prev = v_curr;
	}

	return 0;
}

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

int seminar()
{
	// Create 3 vertices, 1 simplex, 1 mesh
	Cartosphere::Point A(Cartosphere::Image(1, 0, 0));
	Cartosphere::Point B(Cartosphere::Image(0, 1, 0));
	Cartosphere::Point C(Cartosphere::Image(0, 0, 1));
	Cartosphere::Triangle T(A, B, C);
	Cartosphere::TriangularMesh mesh(T);

	// Function
	auto f = [](const Cartosphere::Point& p) -> FLP {
		return 1 - M_2_PI * p.p();
	};

	// Compute centroid/simpsons/3V and iteratively perform mid-point rule
	std::vector<FLP> approximation;
	for (size_t i = 0; i <= 10; ++i)
	{
		approximation.push_back(mesh.integrate(f,
			Cartosphere::TriangularMesh::Quadrature::AreaWeighted,
			Cartosphere::Triangle::Integrator::ThreeVertices));

		if (i < 10)
			mesh.refine();
	}

	// Output results
	for (size_t i = 0; i < approximation.size(); ++i)
	{
		FLP e = std::abs(approximation[i] - (M_PI_2 - 1));
		std::cout << "Refinement " << i << ": "
			<< approximation[i] << " & "
			<< e << " & "
			<< std::log10(e) << '\n';
	}

	return 0;
}

int convergence()
{
	// Files to loop through
	std::vector<std::string> paths = {
		"icosahedron.csm",
		"icosahedron.csm.r1",
		"icosahedron.csm.r2",
		"icosahedron.csm.r3",
		"icosahedron.csm.r4",
		"icosahedron.csm.r5",
	};

	// The function u(x,y) = x y
	auto u = [](const Cartosphere::Point& p) -> FLP
	{
		return pow(p.x(), 2) + pow(p.y(), 2);
	};

	// Laplace-Beltrami of u(x,y,z) = x y
	/*auto f = [](const Cartosphere::Point& p) -> FLP {
		FLP x = p.x(), y = p.y(), z = p.z();
		return x - 2 * x * x * x + y - 2 * x * y - 4 * x * x * y + 6 * x * x * x * y
			- 4 * x * y * y - 2 * y * y * y + 6 * x * y * y * y
			- 2 * x * y * z - 2 * x * z * z - 2 * y * z * z + 6 * x * y * z * z;
	};*/

	auto f = [](const Cartosphere::Point& p) -> FLP {
		FLP x = p.x(), y = p.y(), z = p.z();
		FLP l = 4 - 10 * pow(x, 2) - 10 * pow(y, 2) - (
			x * (x * (2 - 6 * pow(x, 2) - 2 * pow(y, 2)) + y * (-4 * x * y)) +
			y * (x * (-4 * x * y) + y * (2 - 2 * pow(x, 2) - 6 * pow(y, 2))) +
			z * (x * (-4 * x * z) + y * (-4 * y * z) + z * (-2 * pow(x, 2) - 2 * pow(y, 2)))
			);
		return -l;
	};

	int vvv = 4;
	auto ttt = 4;

	// Calculate errors for x = A^{-1}b for each mesh
	std::vector<FLP> L2E;
	for (auto& path : paths)
	{
		FLP error;

		// Construct mesh and linear system
		Cartosphere::TriangularMesh mesh(path);
		Matrix A;
		Vector b;
		mesh.fill(A);
		mesh.fill(b, f);

		// Solve the linear system
		Solver s(A);
		Vector x = s.solve(b);

		// Evaluate the exact solution at the given vertices
		auto vertices = mesh.vertices();
		std::vector<FLP> us;
		std::transform(vertices.cbegin(), vertices.cend(),
			std::back_inserter(us), u);

		// Calculate the ell2-norm of the residual
		FLP norm = 0;
		FLP norm_max = 0;
		std::vector<FLP> e(us.size());
		for (size_t k = 0; k < us.size(); ++k)
		{
			e[k] = us[k] - x[k];
			norm += std::pow(e[k], 2);
			if (abs(e[k]) > norm_max)
			{
				norm_max = abs(e[k]);
			}
			e[k] = pow(e[k], 2);
		}
		error = std::sqrt(norm) / us.size();

		std::cout << "L_inf = " << norm_max << ", L_2 = " << sqrt(mesh.integrate(e)) << "\n";

		// Save error
		L2E.push_back(error);

		// Output colored polyhedral
		std::string output_name;
		{
			std::stringstream sst;
			sst << path << ".sol.obj";
			output_name = sst.str();
		}
		std::vector<FLP> solution;
		solution.reserve(x.size());
		for (int k = 0; k < x.size(); ++k)
		{
			solution.push_back(x[k]);
		}
		mesh.formatPoly(output_name, solution);
	}

	for (size_t k = 0; k < L2E.size(); ++k)
	{
		std::cout << "k=" << k << ", R=" << L2E[k] << "\n";
	}

	return 0;
}

int precompute_weights(const std::string& path)
{
	Cartosphere::TriangularMesh mesh(path);
	// TODO
	return 0;
}

int refine(const std::string& path)
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

int test_obj()
{
	// Let's test the coloring of f(r,theta) = theta!
	// First, we create a mesh
	Cartosphere::TriangularMesh m;
	m.load("icosahedron.csm");

	// Obtain vertices and color linearly
	auto vs = m.vertices();
	std::vector<FLP> xs;

	for (auto& v : vs)
	{
		FLP x = v.x();
		xs.push_back(x);
	}

	// Generate mesh and view in meshlab
	m.format("icosahedron-x-linear.obj", xs);

	return 0;
}

int research_a()
{
	// Set mesh, equation, refinement levels
	std::string name = "icosahedron.csm";
	const int scenario = 0;
	const int refinements = 6;

	// Load initial mesh from path and print statistics
	Cartosphere::TriangularMesh m(name);
	if (!m.isReady())
	{
		// Print error messages
		for (auto& msg : m.getMessages())
		{
			std::cout << msg << "\n";
		}
		std::cout << std::flush;
		return -1;
	}

	auto stats = m.statistics();
	size_t euler = stats.V + stats.F - stats.E;

	std::cout << "Loaded mesh from file: " << name << "\n\n"
		<< "Statistics:\n"
		<< "    Euler: V - E + F = " << stats.V << " - " << stats.E
		<< " + " << stats.F << " = " << euler << "\n"
		<< "    Area ratio: " << stats.areaElementDisparity
		<< " (max " << stats.areaElementMax
		<< ", min " << stats.areaElementMin << ")\n" << std::endl;

	Cartosphere::Function u_inf_func, f_func, g_func;
	if (scenario == 0)
	{
		// The desired steady-state solution u = x^2 + y^2 - 2/3
		// The external term                 f = -Lapl u
		// The initial condition             g = 0
		u_inf_func = [](const Cartosphere::Point& p) -> FLP
		{
			FLP x = p.x(), y = p.y(), z = p.z();
			return pow(x, 2) + pow(y, 2) - FLP(2) / 3;
		};
		f_func = [](const Cartosphere::Point& p) -> FLP {
			FLP x = p.x(), y = p.y(), z = p.z();
			FLP l = 4 - 10 * pow(x, 2) - 10 * pow(y, 2) - (
				x * (x * (2 - 6 * pow(x, 2) - 2 * pow(y, 2)) + y * (-4 * x * y)) +
				y * (x * (-4 * x * y) + y * (2 - 2 * pow(x, 2) - 6 * pow(y, 2))) +
				z * (x * (-4 * x * z) + y * (-4 * y * z) + z * (-2 * pow(x, 2) - 2 * pow(y, 2)))
				);
			return -l;
		};
		g_func = [](const Cartosphere::Point& p) -> FLP {
			FLP x = p.x(), y = p.y(), z = p.z();
			return 0;
		};
	}
	else
	{
		// The desired steady-state solution u = 2
		// The external term                 f = 0
		// The initial condition             g = 2 + z
		u_inf_func = [](const Cartosphere::Point& p) -> FLP
		{
			FLP x = p.x(), y = p.y(), z = p.z();
			return 2;
		};
		f_func = [](const Cartosphere::Point& p) -> FLP {
			FLP x = p.x(), y = p.y(), z = p.z();
			return 0;
		};
		g_func = [](const Cartosphere::Point& p) -> FLP {
			FLP x = p.x(), y = p.y(), z = p.z();
			return 2 + z;
		};
	}

	// Iteratively refine and observe the convergence
	for (int i = 0; i <= refinements; ++i, m.refine())
	{
		// Extract vertices
		const auto vs = m.vertices();

		// Initialize the scalar field at t=0 and infinity
		Vector u_init(vs.size());
		Vector u_inf(vs.size());
		std::transform(vs.begin(), vs.end(), u_init.begin(), g_func);
		std::transform(vs.begin(), vs.end(), u_inf.begin(), u_inf_func);

		// Build the linear system
		Matrix A, M;
		Vector F;
		m.fill(A, M);
		m.fill(F, f_func);

		// Perform time-stepping
		const int time_steps = 200;
		const FLP time_elapsed = 10;
		int iteration = 0;
		FLP indicator = 1;
		Vector u_prev = u_init;
		Vector u_curr;
		for (int step = 0; step < time_steps; ++step, u_prev = u_curr)
		{
			FLP duration = time_elapsed / time_steps;

			Matrix LHS = A + M / duration;
			Vector RHS = F + M * u_prev / duration;

			Solver s(LHS);
			u_curr = s.solve(RHS);
		}

		// Report the L2-error by comparing the linearly weighted approximation
		// to the exact steady-state solution
		std::vector<FLP> u;
		std::transform(u_curr.begin(), u_curr.end(), std::back_inserter(u),
			[](auto& u) -> FLP { return FLP(u); }
		);

		indicator = m.lebesgue(u, u_inf_func);
		std::cout << "R" << i
			<< ": h = " << m.statistics().diameterElementMax <<
			" L2e_" << i << " = " << indicator << "\n";
	}

	return 0;
}

int research_b()
{
	Cartosphere::TriangularMesh m;
	m.load("icosahedron.csm");

	if (!m.isReady())
	{
		std::cout << "The mesh is not ready!\n";
		return -1;
	}

	for (size_t k = 0; k <= 7; ++k, m.refine())
	{
		std::vector<FLP> values(m.vertices().size(), 1);
		FLP integral = m.integrate(values);

		std::cout << "I_" << k << " = " << integral << "\n";
	}

	return 0;
}

int research_c(int l, int m, bool silent)
{
	// Set levels of refinements
	const int levels = 6;

	// Set name of mesh file to load
	std::string name = "icosahedron.csm";

	// Load mesh for steady state solver
	Cartosphere::TriangularMesh mesh;

	mesh.load(name);
	if (!mesh.isReady())
	{
		std::cout << "The mesh is not ready!\n";
		return -1;
	}

	// Print mesh statistics
	auto stats = mesh.statistics();
	size_t euler = stats.V + stats.F - stats.E;

	if (!silent)
	{
		std::cout << "Loaded mesh from file: " << name << "\n"
			<< "Setting up the solver...\n";
	}

	// Set up iterative refinement and the solver
	Cartosphere::SteadyStateSolver solver;
	std::vector<FLP> errors(levels + 1);
	for (int level = 0; level <= levels; ++level, mesh.refine())
	{
		// Output updated mesh information
		stats = mesh.statistics();

		if (!silent)
		{
			std::cout << "Refinement level " << level << "\n"
				<< "Statistics:\n"
				<< "    Euler: V - E + F = " << stats.V << " - " << stats.E
				<< " + " << stats.F << " = " << euler << "\n"
				<< "    Area ratio: " << stats.areaElementDisparity
				<< " (max " << stats.areaElementMax
				<< ", min " << stats.areaElementMin << ")\n"
				<< "    Max diameter: " << stats.diameterElementMax << "\n\n";
		}

		// Update solver
		solver.set(mesh);

		// *******************************************************************
		// Equation                          -LAPL u = f
		// The desired steady-state solution u = Y_l^m
		// The external term                 f = Y_l^m * (l * (l + 1))
		// *******************************************************************
		// When l = 0, u is not zero-averaged
		// When l > 0, u is zero-averged

		Cartosphere::Function u = [l, m](const Cartosphere::Point& p) -> FLP {
			return cartosphere_Y_real(l, m, p.p(), p.a());
		};

		Cartosphere::Function f = [l, m](const Cartosphere::Point& p) -> FLP {
			return l * (l + 1) * cartosphere_Y_real(l, m, p.p(), p.a());
		};

		// Debug the system
		std::string debug_name;
		{
			std::stringstream sst;
			sst << "system_" << l
				<< (m < 0 ? "_m" : "_") << std::abs(m)
				<< "_" << level << ".m";
			debug_name = sst.str();
		}
		// solver.debug(debug_name);

		// Solve the system
		solver.solve(f);
		std::vector<FLP> solution = solver.get();

		// Gauge the error
		FLP error = mesh.lebesgue(solution, u);
		errors[level] = error;

		if (!silent)
		{
			std::cout
				<< "Refine   = " << level << "\n"
				<< "Diameter = " << stats.diameterElementMax << "\n"
				<< "Error    = " << errors[level] << "\n"
				<< std::endl;
		}
		else
		{
			std::cout << "level=" << level
				<< " h=" << stats.diameterElementMax
				<< " e=" << errors[level] << "\n";
		}
	}

	return 0;
}

int research_d()
{

	if (true)
	{
		Cartosphere::Point A(0, 0);
		Cartosphere::Point B(M_PI / 3, 0);
		Cartosphere::Point C(M_PI / 3, M_PI / 3);
		Cartosphere::Triangle ABC(A, B, C);

		FL3 u = ABC.gradient(0);
		FL3 v = ABC.gradient(1);
		FL3 w = ABC.gradient(2);

		std::cout
			<< "A(" << A.x() << ", " << A.y() << ", " << A.z() << ")\n"
			<< "B(" << B.x() << ", " << B.y() << ", " << B.z() << ")\n"
			<< "C(" << C.x() << ", " << C.y() << ", " << C.z() << ")\n"
			<< "u<" << u.x << ", " << u.y << ", " << u.z << ">\n"
			<< "v<" << v.x << ", " << v.y << ", " << v.z << ">\n"
			<< "w<" << w.x << ", " << w.y << ", " << w.z << ">\n";

		auto D = Cartosphere::Arc(B, C).midpoint();

		std::cout
			<< "D(" << D.x() << ", " << D.y() << ", " << D.z() << ")\n";

		FLP d = Cartosphere::Arc(D, A).length();

		std::cout << "d = " << d << "\n";

		FL3 t = Cartosphere::Arc(D, A).tangent(0);

		std::cout
			<< "t<" << t.x << ", " << t.y << ", " << t.z << ">\n";

		FL3 g = ABC.gradient(0);

		std::cout
			<< "g<" << g.x << ", " << g.y << ", " << g.z << ">\n";

		FL3 T = Cartosphere::transport(A, D, g);

		std::cout
			<< "T<" << T.x << ", " << T.y << ", " << T.z << ">\n";

		std::cout
			<< "2T = <" << (2 * T).x << ", " << (2 * T).y << ", " << (2 * T).z << ">\n";
	}

	{
		Cartosphere::Point B(0, 0);
		Cartosphere::Point C(M_PI / 3, 0);
		Cartosphere::Point A(M_PI / 3, M_PI / 3);
		Cartosphere::Triangle ABC(A, B, C);

		auto D = Cartosphere::Arc(C, A).midpoint();

		FL3 g = ABC.gradient(1);

		FL3 T = Cartosphere::transport(B, D, g);

		std::cout
			<< "T<" << T.x << ", " << T.y << ", " << T.z << ">\n";
	}

	{
		Cartosphere::Point C(0, 0);
		Cartosphere::Point A(M_PI / 3, 0);
		Cartosphere::Point B(M_PI / 3, M_PI / 3);
		Cartosphere::Triangle ABC(A, B, C);

		auto D = Cartosphere::Arc(A, B).midpoint();

		FL3 g = ABC.gradient(2);

		FL3 T = Cartosphere::transport(C, D, g);

		std::cout
			<< "T<" << T.x << ", " << T.y << ", " << T.z << ">\n";
	}


	{
		Cartosphere::Point A(0, 0);
		Cartosphere::Point B(M_PI / 3, 0);
		Cartosphere::Point C(M_PI / 3, M_PI / 3);
		Cartosphere::Triangle ABC(A, B, C);

		std::cout
			<< ABC.contains(A)
			<< ABC.contains(B)
			<< ABC.contains(C)
			<< ABC.contains(Cartosphere::Point(M_PI / 6, M_PI / 6)) << "\n";
	}

	TriangularMesh m("icosahedron.csm.5.csm");
	// m.format("icosahedron.csm.5.csm.obj");

	std::vector<Point> v = m.vertices();
	std::vector<FLP> a;
	std::vector<FLP> t(101, 0);
	for (int k = 0; k < t.size(); ++k)
	{
		t[k] = 0.1 * k;
	}
	
	std::vector<Point> p;
	for (int k = -180; k < 180; ++k)
	{
		p.emplace_back(M_PI_2+FLP(1e-6), deg2rad(k));
	}
	
	TimeDependentSolver s;
	s.set(m);
	s.set([](const Point& x) -> FLP { return 0; });
	s.initialize([](const Point& x) -> FLP { return 2 + x.z(); });
	for (size_t k = 0; k + 1 < t.size(); ++k)
	{
		FLP duration = t[k + 1] - t[k];
		s.advance(duration);
		std::vector<FL3> u = s.velocity(p);
		std::vector<Point> q = p;
		for (size_t i = 0; i < p.size(); ++i)
		{
			p[i].move(duration * u[i]);
		}
		// std::cout << "k = " << k << '\n';

		// std::cout << p[1].x() << " " << p[1].y() << " " << p[1].z() << ";\n";

	}


	std::vector<FLP> ps;

	for (size_t i = 0; i < p.size(); ++i)
	{
		ps.push_back(p[i].p());
		// std::cout << "p[" << i << "] = (" << rad2deg(p[i].p()) << ", " << rad2deg(p[i].a()) << ")\n";
	}

	double sum = std::accumulate(ps.begin(), ps.end(), 0.0);
	double mean = sum / ps.size();

	double sq_sum = std::inner_product(ps.begin(), ps.end(), ps.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / ps.size() - mean * mean);

	std::cout << rad2deg(mean) << " +/- " << rad2deg(stdev) << "\n";

	return 0;
}

int research_f()
{
	TriangularMesh mesh;
	mesh.load("icosahedron.csm");
	if (!mesh.isReady())
	{
		std::cout << "Mesh not ready.\n" << std::endl;
		auto ms = mesh.getMessages();
		for (auto m : ms)
		{
			std::cout << m << "\n";
		}
		return -1;
	}
	std::cout << "Mesh is ready!\n" << std::endl;

	const FLP time_initial = 1e-4;
	const FLP time_ratio = 1.01;
	const FLP dist_tolerance = 1e-7;
	const FLP time_max = 50;

	Cartosphere::TimeDependentSolver solver;
	size_t levels = 6;
	for (int level = 0; level <= levels; ++level, mesh.refine())
	{
		// Output updated mesh information
		auto stats = mesh.statistics();
		size_t euler = stats.V + stats.F - stats.E;

		std::cout << "Refinement level " << level << "\n"
			<< "Statistics:\n"
			<< "    Euler: V - E + F = " << stats.V << " - " << stats.E
			<< " + " << stats.F << " = " << euler << "\n"
			<< "    Area ratio: " << stats.areaElementDisparity
			<< " (max " << stats.areaElementMax
			<< ", min " << stats.areaElementMin << ")\n"
			<< "    Max diameter: " << stats.diameterElementMax << "\n\n";
		
		std::vector<Point> p;
		for (int k = -180; k < 180; ++k)
		{
			p.emplace_back(M_PI_2, deg2rad(k));
		}
		std::vector<Point> targets;
		for (int k = -180; k < 180; ++k)
		{
			targets.emplace_back(acos(-0.25), deg2rad(k));
		}

		// Update solver
		TimeDependentSolver solver;
		solver.set(mesh);
		solver.set([](const Point& x) -> FLP { return 0; });
		solver.initialize([](const Point& x) -> FLP { return 2 + x.z(); });

		size_t step = 0;
		size_t cumulative = 0;
		FLP duration = time_initial;
		FLP dist_max = std::numeric_limits<FLP>::max();
		FLP vertex_change = std::numeric_limits<FLP>::max();
		solver.advance(duration/2);
		while (cumulative < time_max && vertex_change > 1e-6)
		{
			std::vector<FL3> u = solver.velocity(p);
			std::vector<Point> q = p;
			dist_max = 0;
			for (size_t i = 0; i < p.size(); ++i)
			{
				auto q = p[i];
				p[i].move(duration * u[i]);
				FLP dist = u[i].norm2();
				if (dist > dist_max)
				{
					dist_max = dist;
				}
			}

			std::cout << "Iteration " << step
				<< " time " << duration
				<< " disp " << dist_max 
				<< " " << (acos(-0.25) - p[0].p()) << "\n";

			++step;
			cumulative += duration;
			duration *= time_ratio;
			vertex_change = solver.advance(duration);

			// std::cout << vertex_change << '\n';
		}

		std::vector<FLP> ps;
	
		for (size_t i = 0; i < p.size(); ++i)
		{
			ps.push_back(distance(p[i], targets[i]));
			// std::cout << "p[" << i << "] = (" << rad2deg(p[i].p()) << ", " << rad2deg(p[i].a()) << ")\n";
		}
	
		double sum = std::accumulate(ps.begin(), ps.end(), 0.0);
		double mean = sum / ps.size();
		double sq_sum = std::inner_product(ps.begin(), ps.end(), ps.begin(), 0.0);
		double stddev = std::sqrt(sq_sum / ps.size() - mean * mean);
	
		std::cout << mean << " +/- " << stddev << "\n";
	}


	return 0;
}

int research_g(const std::string &folder)
{
	std::cout << "Initializing from directory " << folder << "\n";

	Cartosphere::ShapeFile file;

	std::string message;
	if (!file.open(folder, message))
	{
		std::cout << "Error: " << message << "\n";
		return -1;
	}

	std::cout << "Initialization success!\n";
	
	std::cout << "Shapes loaded: " << file.count() << "\n";

	return 0;
}

int
benchmark()
{
	std::cout << "HERE WE GO AGAIN\n";

	Cartosphere::SpectralSolver s;
	s.parse("population.csm");

	std::cout << s.inputSummary() << std::endl;
	s.execute();

	std::cout << s.outputSummary() << std::endl;

	std::cout << "EVERYTHING PASSED!\n";

	return 0;
}
