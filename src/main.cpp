
#include <argparse/argparse.hpp>
using argparse::ArgumentParser;

#include <iostream>

#include "cartosphere/research.hpp"

int
runDemo(const std::string&, const std::vector<std::string>&);

/* **************************** MAIN ENTRY POINT **************************** */
int
main(int argc, char* argv[])
{
	ArgumentParser args("cartosphere", "0.1.0-dev");

	args.add_argument("--demo")
		.default_value(std::vector{ std::string{"list"}})
		.help("run demo scenarios")
		.nargs(argparse::nargs_pattern::any);

	try
	{
		args.parse_args(argc, argv);
	}
	catch (const std::runtime_error& err)
	{
		std::cerr << err.what() << "\n";
		std::cerr << args;
		std::exit(1);
	}

	auto isDemo = args.is_used("--demo");
	
	if (isDemo)
	{
		auto demoArgs = args.get<std::vector<std::string>>("--demo");
		auto demoName = demoArgs.front();
		demoArgs.erase(demoArgs.begin());

		return runDemo(demoName, demoArgs);
	}

	return 0;
}

int
runDemo(const std::string& name, const std::vector<std::string>& args)
{
	if (name == "list")
	{
		std::cout << "Available demos:\n"
			<< "demo\n"
			<< "diffusion\n"
			<< "seminar\n"
			<< "quadrature\n"
			<< "testobj\n"
			<< "benchmark\n"
			<< "precompute\n"
			<< "refine\n"
			<< "A\n"
			<< "B\n"
			<< "C\n"
			<< "CC\n"
			<< "D\n"
			<< "F\n"
			<< "G\n";
		return 0;
	}

	if (name == "demo")
		return demo();
	
	if (name == "diffusion")
		return demo_diffusion();

	if (name == "seminar")
		return seminar();

	if (name == "quadrature")
		return demo_quadrature();

	if (name == "testobj")
		return test_obj();

	if (name == "benchmark")
		return benchmark();

	if (name == "precompute")
	{
		if (args.size() != 1)
			throw std::runtime_error("This demo needs 1 argument.");

		auto path = args[0];
		return precompute_weights(path);
	}

	if (name == "refine")
	{
		if (args.size() != 1)
		{
			std::cerr << "Needs 1 demo argument.";
			return 1;
		}

		auto path = args[0];
		return refine(path);
	}

	// RESEARCH CODE
	if (name == "A")
		return research_a();

	if (name == "B")
		return research_b();

	if (name == "C")
	{
		if (args.size() != 2)
		{
			std::cerr << "Needs 2 demo arguments.";
			return 1;
		}

		auto l = std::stoi(args[0]);
		auto m = std::stoi(args[1]);
		return research_c(l, m);
	}

	if (name == "CC")
	{
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

	if (name == "D")
		return research_d();

	if (name == "F")
		return research_f();

	if (name == "G")
	{
		if (args.size() != 1)
		{
			std::cerr << "Needs 1 demo argument.\n";
			return 1;
		}

		auto folder = args[0];
		return research_g(folder);
	}

	std::cerr << "Unknown demo name\n";
	return 0;
}
