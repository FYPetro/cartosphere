
#include "cartosphere/research.hpp"

#include <iostream>

/* **************************** MAIN ENTRY POINT **************************** */
int main(int argc, char** argv)
{
	// The parsing of arguments
	if (argc < 2)
	{
		return demo_diffusion();
	}

	// Switch based on the first keyword
	std::string item(argv[1]);
	if (item == "demo")
	{
		return demo();
	}
	else if (item == "seminar")
	{
		return seminar();
	}
	else if (item == "quadrature")
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
	else if (item == "convergence")
	{
		return convergence();
	}
	else if (item == "testobj")
	{
		return test_obj();
	}
	else if (item == "research")
	{
		if (argc < 3)
		{
			std::cout << "Specify research option with res a/b/c/etc\n";
			return 0;
		}

		std::string option(argv[2]);
		if (option == "a")
		{
			// Find the rate of convergence of the diffusion problem
			return research_a();
		}
		else if (option == "b")
		{
			// Check aall error gauging functions
			return research_b();
		}
		else if (option == "c")
		{
			// Test steady-state solver for FEM
			if (argc < 5)
			{
				std::cout << "Specify degree and order\n";
				return 0;
			}
			return research_c(std::stoi(argv[3]), std::stoi(argv[4]));
		}
		else if (option == "cc")
		{
			// Test all 0 < l <= 1
			for (int l = 1; l <= 3; ++l)
			{
				std::cout << "Y_" << l << "^" << 0 << "\n\n";
				research_c(l, 0, true);

				for (int m = 1; m <= l; ++m)
				{
					std::cout << "Y_" << l << "^" << m << "\n\n";
					research_c(l, m, true);

					std::cout << "Y_" << l << "^" << -m << "\n\n";
					research_c(l, -m, true);
				}
			}
			return 0;
		}
		else if (option == "d")
		{
			return research_d();
		}
		else if (option == "f")
		{
			return research_f();
		}
	}
	else if (item == "benchmark")
	{
		return benchmark();
	}

	return 0;
}
