
#include <argparse/argparse.hpp>
using argparse::ArgumentParser;

#include <iostream>

#include "cartosphere/shapefile.hpp"

#include "cartosphere/research.hpp"

int
runDemo(const std::string&, const std::vector<std::string>&);

/* **************************** MAIN ENTRY POINT **************************** */
int
main(int argc, char* argv[])
{
	ArgumentParser args("cartosphere", "0.1.0-dev");

	args.add_argument("mode")
		.help("specify cartosphere operation mode")
		.default_value(std::string{ "demo" })
#if __APPLE__
		.metavar("MODE")
#endif
		;
	// cartosphere viz -i cb_2021_us_state_500k -o cb_2021_us_state_500k.m

	// Demo overrides all other options
	args.add_argument("-s", "--scene")
		.help("specify demo scene")
		.nargs(argparse::nargs_pattern::at_least_one)
		.default_value(std::vector{ std::string{"list"}})
#if __APPLE__
		.metavar("SCENE [args...]")
#endif
		;

	// If not a demo, an input must be specified
	args.add_argument("-i", "--input")
		.help("path of input file/folder")
		.nargs(1)
#if __APPLE__
		.metavar("INPUT")
#endif
		;

	args.add_argument("--input-format")
		.help("format of input")
		.default_value(std::string{ "shapefile" })
#if __APPLE__
		.metavar("INFMT")
#endif
		;

	args.add_argument("-o", "--output")
		.help("path to output file/folder")
		.nargs(1)
#if __APPLE__
		.metavar("OUTPUT")
#endif
		;

	args.add_argument("--output-format")
		.help("format of output")
		.default_value(std::string{ "matlab" })
#if __APPLE__
		.metavar("OUTFMT")
#endif
		;

	args.add_argument("-m", "--mesh")
		.help("specify input .csm file as background mesh")
		.nargs(1)
#if __APPLE__
		.metavar("CSMFILE")
#endif
		;

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

	auto mode = args.get<std::string>("mode");

	// Priotize demo over all other operating modes
	if (mode == "demo")
	{
		auto demoArgs = args.get<std::vector<std::string>>("--scene");
		auto demoName = demoArgs.front();
		demoArgs.erase(demoArgs.begin());

		return runDemo(demoName, demoArgs);
	}

	// Visualize a file
	if (mode == "viz")
	{
		auto inputPath = args.get<std::string>("--input");
		auto inputFormat = args.get<std::string>("--input-format");
		std::cout << "Input path: " << inputPath
			<< " (format: " << inputFormat << ")\n";

		auto outputPath = args.get<std::string>("--output");
		auto outputFormat = args.get<std::string>("--output-format");
		std::cout << "Output path: " << outputPath
			<< " (format: " << outputFormat << ")\n";

		// Visualize a shapefile
		if (inputFormat == "shapefile")
		{
			std::cout << "Initializing shapefile from " << inputPath << "...\n";

			Cartosphere::ShapeFile shapefile;
			std::string message;
			if (!shapefile.open(inputPath, message))
			{
				std::cerr << "Error: " << message << "\n";
				std::exit(1);
			}
			std::cout << "Shapes loaded: " << shapefile.count() << "\n";

			if (outputFormat == "matlab")
			{
				std::cout << "Viz-ing shapefile using matlab...\n";
				std::exit(0);
			}

			std::cerr << "Unhandled Output Format: " << outputFormat << "\n";
			std::exit(1);
		}

		std::cerr << "Unhandled Input Format: " << inputFormat << "\n";
		std::exit(1);
	}

	// Command not recognized
	std::cerr << "Unrecognized mode " << mode << "\n";
	std::exit(1);

	return 0;
}

int
runDemo(const std::string& name, const std::vector<std::string>& args)
{
	if (name == "list")
	{
		std::cout << "Available demo SCENARIOs:\n"
			<< "\tdemo               [---]\n"
			<< "\tdiffusion          [---]\n"
			<< "\tseminar            [---]\n"
			<< "\tquadrature         [---]\n"
			<< "\ttestobj            [---]\n"
			<< "\tbenchmark          [---]\n"
			<< "\tprecompute         [---]\n"
			<< "\trefine LEVEL       [---]\n"
			<< "\tA                  [Research A]\n"
			<< "\tB                  [Research B]\n"
			<< "\tC L M              [Research C]\n"
			<< "\tCC                 [Research CC]\n"
			<< "\tD                  [Research D]\n"
			<< "\tF                  [Research F]\n"
			<< "\tG SHAPEFILE        [Research G]\n"
			<< "Usage: cartosphere demo --scene [SCENARIO]\n";
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
