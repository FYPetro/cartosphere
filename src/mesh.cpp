
#include "cartosphere/predicates.hpp"
#include "cartosphere/mesh.hpp"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace Cartosphere::Predicates;

Cartosphere::Image
Cartosphere::Preimage::toImage() const
{
	Image image;

	FLP projection = sin(p);

	image.x = projection * cos(a);
	image.y = projection * sin(a);
	image.z = cos(p);

	return image;
}

Cartosphere::Preimage
Cartosphere::Image::toPreimage() const
{
	Preimage preimage;

	preimage.p = acos(z);

	if (abs(x) + abs(y) > EPS)
	{
		preimage.a = atan2(x, y);
	}
	else
	{
		preimage.a = 0;
	}

	return preimage;
}

FLP
Cartosphere::Triangle::area() const
{
	return 4;
}

FLP
Cartosphere::Triangle::areaEuclidean() const
{
	// Form vector AB
	FL3 AB = B.image() - A.image();
	FL3 AC = C.image() - A.image();
	FL3 product = cross(AB, AC);
	// Form vector AC
	return 0.5 * product.norm2();
}

void
Cartosphere::TriangularMesh::clear()
{
	// Reset state variables
	listofPoints.clear();
	listofEdges.clear();
	listofTriangles.clear();
	listofMessages.clear();

	// Reset state flags
	flagofLoading = false;
}

bool
Cartosphere::TriangularMesh::load(std::string const &path)
{
	if (flagofLoading) clear();

	// Opens given path
	std::ifstream ifs(path);
	if (!ifs.is_open())
	{
		flagofLoading = false;
		std::string message;
		{
			std::stringstream sst;
			sst << "Coud not load path " << path;
			message = sst.str();
		}
		listofMessages.push_back(message);
		return isLoaded();
	}

	// Temporary state variables
	unsigned lineNumber = 0;
	unsigned lineParsed = 0;
	std::vector<unsigned> specs;
	FLP coords[4];
	// Temporary state flags
	unsigned format = 0;

	// Start parsing file line by line
	std::string line;
	std::istringstream iss;
	while (std::getline(ifs, line))
	{
		// Increment line number and skip empty lines
		++lineNumber;
		if (!line.length()) continue;
		iss.str(line);

		// Skip comment lines and strip inline comments
		if (line[0] == '#') continue;
		unsigned indexofPound = line.find_first_of('#');
		if (indexofPound != std::string::npos)
		{
			line.erase(indexofPound);
		}

		// Parse sanitized inputs
		if (lineParsed == 0)
		{
			// Section 0: Size specifications
			unsigned size;
			
			iss >> size;
			if (iss.fail() || size == 0)
			{
				std::string message;
				{
					std::stringstream sst;
					sst << "Error in Line " << lineNumber
						<< ": Number of points is missing or zero";
					message = sst.str();
				}
				listofMessages.push_back(message);
				return isLoaded();
			}
			specs.push_back(size);

			iss >> size;
			if (iss.fail() || size == 0)
			{
				std::string message;
				{
					std::stringstream sst;
					sst << "Error in Line " << lineNumber
						<< ": Number of edges is missing or zero";
					message = sst.str();
				}
				listofMessages.push_back(message);
				return isLoaded();
			}
			specs.push_back(size);

			iss >> size;
			if (iss.fail() || size == 0)
			{
				std::string message;
				{
					std::stringstream sst;
					sst << "Error in Line " << lineNumber
						<< ": Number of triangles is missing or zero";
					message = sst.str();
				}
				listofMessages.push_back(message);
				return isLoaded();
			}
			specs.push_back(size);

			iss >> format;
			if (iss.fail())
			{
				format = 0;
			}
			else if (format != 0)
			{
				std::string message;
				{
					std::stringstream sst;
					sst << "Error in Line " << lineNumber
						<< ": File format ID " << format << " is unrecognized";
					message = sst.str();
				}
				listofMessages.push_back(message);
				return isLoaded();
			}
			specs.push_back(format);
		}
		else if (lineParsed <= specs[0])
		{
			// Section 1: List of points
			for (unsigned pos = 0; pos < 2; ++pos)
			{
				iss >> coords[pos];
				if (iss.fail())
				{
					std::string message;
					{
						std::stringstream sst;
						sst << "Error in Line " << lineNumber
							<< ": Missing coordinate " << pos;
						message = sst.str();
					}
					listofMessages.push_back(message);
					return isLoaded();
				}
			}

			iss >> coords[2];
			if (iss.fail())
			{
				// Spherical coordinates in degrees
				Preimage pi(deg2rad(coords[0]), deg2rad(coords[1]));
				listofPoints.emplace_back(pi);
			}
			else
			{
				// Cartesian coordinates in degrees
				Image im(coords[0], coords[1], coords[2]);
				listofPoints.emplace_back(im);
			}

			iss >> coords[3];
			if (!iss.fail())
			{
				std::string message;
				{
					std::stringstream sst;
					sst << "Warning in Line " << lineNumber
						<< ": Extra arguments are dropped";
					message = sst.str();
				}
				listofMessages.push_back(message);
			}
		}
		else if (lineParsed <= specs[0] + specs[1])
		{
			// Section 2: List of edges
			std::pair<unsigned, unsigned> edge;
			iss >> edge.first;
			iss >> edge.second;
			if (iss.fail())
			{
				std::string message;
				{
					std::stringstream sst;
					sst << "Error in Line " << lineNumber
						<< ": Edge specification missing point(s) ";
					message = sst.str();
				}
				listofMessages.push_back(message);
				return isLoaded();
			}
			listofEdges.push_back(edge);

			iss >> edge.first;
			if (!iss.fail())
			{
				std::string message;
				{
					std::stringstream sst;
					sst << "Warning in Line " << lineNumber
						<< ": Extra arguments are dropped";
					message = sst.str();
				}
				listofMessages.push_back(message);
			}
		}
		else if (lineParsed <= specs[0] + specs[1] + specs[2])
		{
			// Section 3: List of triangles
			unsigned indices[3];
			char orientation[3] = { '\0', '\0', '\0' };

			for (unsigned pos = 0; pos < 3; ++pos)
			{
				iss >> orientation[pos];
				if (iss.fail() || orientation[0] == '\0')
				{
					std::string message;
					{
						std::stringstream sst;
						sst << "Error in Line " << lineNumber
							<< ": Argument " << pos << " is missing an orientation";
						message = sst.str();
					}
					listofMessages.push_back(message);
				}
				iss >> indices[pos];
				if (iss.fail() || orientation[0] == '\0')
				{
					std::string message;
					{
						std::stringstream sst;
						sst << "Error in Line " << lineNumber
							<< ": Argument " << pos << " is not formatted correctly";
						message = sst.str();
					}
					listofMessages.push_back(message);
				}
			}

			// Construct and push triangle
			Triangle triangle;
			if (orientation[0] == '+')
				triangle.A = listofPoints[listofEdges[indices[0]].first];
			else
				triangle.A = listofPoints[listofEdges[indices[0]].second];
			if (orientation[1] == '+')
				triangle.B = listofPoints[listofEdges[indices[1]].first];
			else
				triangle.B = listofPoints[listofEdges[indices[1]].second];
			if (orientation[2] == '+')
				triangle.C = listofPoints[listofEdges[indices[2]].first];
			else
				triangle.C = listofPoints[listofEdges[indices[2]].second];
			listofTriangles.push_back(triangle);

			iss >> orientation[0];
			if (!iss.fail())
			{
				std::string message;
				{
					std::stringstream sst;
					sst << "Warning in Line " << lineNumber
						<< ": Extra arguments are dropped";
					message = sst.str();
				}
				listofMessages.push_back(message);
			}
		}

		++lineParsed;
		iss.clear();
	}

	flagofLoading = true;
	return isLoaded();
}

void
Cartosphere::TriangularMesh::reportAreas()
{
	FLP area = 0;
	FLP totalArea = 0;

	for (auto &triangle : listofTriangles)
	{
		area = triangle.areaEuclidean();
		totalArea += area;
		std::string message;
		{
			auto &A = triangle.A.image();
			auto &B = triangle.B.image();
			auto &C = triangle.C.image();
			std::stringstream sst;
			sst << "Area: " << area
				<< " | A(" << A.x << ", " << A.y << ", " << A.z << ") "
				<< " B(" << B.x << ", " << B.y << ", " << B.z << ") "
				<< " C(" << C.x << ", " << C.y << ", " << C.z << ") ";
			message = sst.str();
		}
		listofMessages.push_back(message);
	}

	std::string message;
	{
		std::stringstream sst;
		sst << "Total Area: " << totalArea;
		message = sst.str();
	}
	listofMessages.push_back(message);
}
