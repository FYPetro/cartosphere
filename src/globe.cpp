
#include "cartosphere/functions.hpp"

#include "cartosphere/globe.hpp"
using Cartosphere::Globe;
using Cartosphere::Polygon;

Globe::Globe(const string& path)
{
	vector<string> lines;
	ifstream ifs(path);
	if (!ifs.is_open())
		return;

	// Fill in the following items
	vector<string> names;
	vector<double> values;
	vector<Polygon> polygons;

	// Read then parse the file line by line
	string line;
	string block;
	string token;
	istringstream sst;
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
				string name;
				sst >> name;
				names.push_back(name);
			}
			else if (token == "VALUE")
			{
				double value;
				sst >> value;
				values.push_back(value);
			}
			else if (token == "PAR")
			{
				double polar, azimuth;
				sst >> polar >> azimuth;
				polygon.emplace_back(polar, azimuth);
			}
			else if (token == "PAD")
			{
				double polar, azimuth;
				sst >> polar >> azimuth;
				polar = cs_deg2rad(polar);
				azimuth = cs_deg2rad(polar);
				polygon.emplace_back(polar, azimuth);
			}
			else if (token == "XYZ")
			{
				double x, y, z;
				sst >> x >> y >> z;
				polygon.emplace_back(x, y, z);
			}
		}
	}

	// Create the mesh
	_m.load("icosahedron.5.csm");

	// Find if vertex is in polygon.
	// vector<Point> v = _m.vertices();
	// vector<double> init(0, v.size());
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

void Globe::format(const string& path)
{
	// Nothing to see here.
}
