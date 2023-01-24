
#include "cartosphere/globe.hpp"
#include "cartosphere/utility.hpp"

#include <vector>
#include <fstream>

using Cartosphere::Globe;
using Cartosphere::Polygon;

Globe::Globe(const std::string& path)
{
	std::vector<std::string> lines;
	std::ifstream ifs(path);
	if (!ifs.is_open())
		return;

	// Fill in the following items
	std::vector<std::string> names;
	std::vector<FLP> values;
	std::vector<Polygon> polygons;

	// Read then parse the file line by line
	std::string line;
	std::string block;
	std::string token;
	std::istringstream sst;
	Polygon polygon;
	while (!ifs.eof())
	{
		std::getline(ifs, line);
		lines.push_back(line);
	}
	for (int k = 0; k < lines.size(); ++k)
	{
		line = trim_copy(lines[k]);
		if (line.empty())
		{
			continue;
		}
		// Exact keyword
		token = line.substr(0, line.find(' '));
		if (block.empty())
		{
			block = token;
		}
		else if (token == "END")
		{
			if (block == "REGION")
			{
				// Save polygon
				polygons.push_back(polygon);
				polygon.clear();
			}
			block.clear();
		}
		else if (block == "REGION")
		{
			sst.str(line.substr(token.size()));
			if (token == "NAME")
			{
				std::string name;
				sst >> name;
				names.push_back(name);
			}
			else if (token == "VALUE")
			{
				FLP value;
				sst >> value;
				values.push_back(value);
			}
			else if (token == "PAR")
			{
				FLP polar, azimuth;
				sst >> polar >> azimuth;
				polygon.emplace_back(polar, azimuth);
			}
			else if (token == "PAD")
			{
				FLP polar, azimuth;
				sst >> polar >> azimuth;
				polar = deg2rad(polar);
				azimuth = deg2rad(polar);
				polygon.emplace_back(polar, azimuth);
			}
			else if (token == "XYZ")
			{
				FLP x, y, z;
				sst >> x >> y >> z;
				polygon.emplace_back(x, y, z);
			}
		}
	}

	// Create the mesh
	_m.load("icosahedron.5.csm");

	// Find if vertex is in polygon.
	// std::vector<Point> v = _m.vertices();
	// std::vector<FLP> init(0, v.size());
	// for (int i = 0; i < v.size(); ++i)
	// {
	// 	for (int j = 0; j < polygons.size(); ++j)
	// 	{
	// 		if (polygons.interior(v[i]))
	// 		{
	// 			init[i] = values[j];
	// 		}
	// 	}
	// }
}

void Globe::run()
{
	// Also nothing to see here.
}

void Globe::format(const std::string& path)
{
	// Nothing to see here.
}
