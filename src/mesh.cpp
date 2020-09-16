
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
		preimage.a = atan2(y, x);
	}
	else
	{
		preimage.a = 0;
	}

	return preimage;
}

Cartosphere::Point
Cartosphere::midpoint(Point const &a, Point const &b)
{
	if (a.isAntipodalTo(b))
	{
		return Point();
	}

	FL3 middleofChord = (FL3)a.image() + (FL3)b.image();
	middleofChord /= 2;
	middleofChord.normalize();
	
	return Point(Image(middleofChord));
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
	listofSimplices.clear();
	listofMessages.clear();

	// Reset state flags
	flagofLoading = false;
	flagofParsing = false;
}

FLP
Cartosphere::TriangularMesh::areaEuclidean() const
{
	FLP area = (FLP)0.0;
	for (auto &simplex : listofSimplices)
	{
		area += simplex.areaEuclidean();
	}
	return area;
}

bool
Cartosphere::TriangularMesh::load(std::string const &path)
{
	if (flagofLoading || flagofParsing) clear();

	// Opens given path
	std::ifstream ifs(path);
	if (!ifs.is_open())
	{
		std::string message;
		{
			std::stringstream sst;
			sst << "Coud not load path " << path;
			message = sst.str();
		}
		listofMessages.push_back(message);
		return isReady();
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
		auto indexofPound = line.find_first_of('#');
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
				return isReady();
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
				return isReady();
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
				return isReady();
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
				return isReady();
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
					return isReady();
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
			PointIndexPair edge;
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
				return isReady();
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

			listofTriangles.emplace_back(
				std::make_pair(indices[0], orientation[0] == '+'),
				std::make_pair(indices[1], orientation[1] == '+'),
				std::make_pair(indices[2], orientation[2] == '+'));

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
	fillSimplices();
	return isReady();
}

bool
Cartosphere::TriangularMesh::save(std::string const &path) const
{
	std::ofstream ofs(path);

	if (!ofs.is_open())
	{
		return false;
	}

	ofs << "#Cartosphere Mesh Fomat\n"
		<< "#V #E #F\n"
		<< listofPoints.size() << " "
		<< listofEdges.size() << " "
		<< listofTriangles.size() << "\n\n"
		<< "#V List\n";

	for (auto &point : listofPoints)
	{
		ofs << rad2deg(point.p()) << " "
			<< rad2deg(point.a()) << "\n";
	}

	ofs << "\n"
		<< "#E List\n";

	for (auto &edge : listofEdges)
	{
		ofs << edge.first << " " << edge.second << "\n";
	}

	ofs << "\n"
		<< "#F List\n";

	for (auto &triangle : listofTriangles)
	{
		auto &a = std::get<0>(triangle);
		auto &b = std::get<1>(triangle);
		auto &c = std::get<2>(triangle);

		ofs << (a.second ? "+" : "-") << a.first << " "
			<< (b.second ? "+" : "-") << b.first << " "
			<< (c.second ? "+" : "-") << c.first << "\n";
	}

	return true;
}

void
Cartosphere::TriangularMesh::refine()
{
	// Reserve space for new midpoints, edges (all), and triangles (all)
	// - listofPoints will accept new midpoints
	// - listofEdges will be cleared
	// - listofTriangles will also be cleared
	listofPoints.reserve(listofPoints.size() + listofEdges.size());
	auto edges = std::move(listofEdges);
	listofEdges.reserve(2 * listofEdges.size() + 3 * listofTriangles.size());
	auto triangles = std::move(listofTriangles);
	listofTriangles.reserve(4 * listofTriangles.size());

	// Each old edge will produce one new midpoint
	for (auto &pair : edges)
	{
		Point const &A = listofPoints[pair.first];
		Point const &B = listofPoints[pair.second];

		unsigned midpointIndex = listofPoints.size();
		listofPoints.push_back(midpoint(A, B));

		// Each new midpoint will bisect the old edge into two new edges
		listofEdges.emplace_back(std::make_pair(pair.first, midpointIndex));
		listofEdges.emplace_back(std::make_pair(midpointIndex, pair.second));
	}

	// Assemble new list of triangles
	unsigned myMidpoints[3];
	std::vector<EdgeIndex> myEdges;
	std::vector<EdgeIndexTriplet> myTriangles;
	myEdges.reserve(12);
	myTriangles.reserve(4);
	for (auto &triangle : triangles)
	{
		auto &a = std::get<0>(triangle);
		auto &b = std::get<1>(triangle);
		auto &c = std::get<2>(triangle);
		// Push 6 constructed edges and retrieve midpoints
		myEdges.emplace_back(2 * a.first, a.second);
		myEdges.emplace_back(2 * a.first + 1, a.second);
		if (!a.second)
			std::swap(myEdges[0], myEdges[1]);

		myEdges.emplace_back(2 * b.first, b.second);
		myEdges.emplace_back(2 * b.first + 1, b.second);
		if (!b.second)
			std::swap(myEdges[2], myEdges[3]);

		myEdges.emplace_back(2 * c.first, c.second);
		myEdges.emplace_back(2 * c.first + 1, c.second);
		if (!c.second)
			std::swap(myEdges[4], myEdges[5]);

		// Retrieve the midpoints
		if (myEdges[0].second)
			myMidpoints[0] = listofEdges[myEdges[0].first].second;
		else
			myMidpoints[0] = listofEdges[myEdges[0].first].first;

		if (myEdges[2].second)
			myMidpoints[1] = listofEdges[myEdges[2].first].second;
		else
			myMidpoints[1] = listofEdges[myEdges[2].first].first;

		if (myEdges[4].second)
			myMidpoints[2] = listofEdges[myEdges[4].first].second;
		else
			myMidpoints[2] = listofEdges[myEdges[4].first].first;

		// Push six new edges from the three new midpoints per triangle
		myEdges.emplace_back(listofEdges.size(), true);
		myEdges.emplace_back(listofEdges.size(), false);
		listofEdges.push_back(std::make_pair(myMidpoints[0], myMidpoints[2]));

		myEdges.emplace_back(listofEdges.size(), true);
		myEdges.emplace_back(listofEdges.size(), false);
		listofEdges.push_back(std::make_pair(myMidpoints[1], myMidpoints[0]));

		myEdges.emplace_back(listofEdges.size(), true);
		myEdges.emplace_back(listofEdges.size(), false);
		listofEdges.push_back(std::make_pair(myMidpoints[2], myMidpoints[1]));

		// Construct four new triangles per old triangle
		listofTriangles.emplace_back(myEdges[0], myEdges[6], myEdges[5]);
		listofTriangles.emplace_back(myEdges[1], myEdges[2], myEdges[8]);
		listofTriangles.emplace_back(myEdges[10], myEdges[3], myEdges[4]);
		listofTriangles.emplace_back(myEdges[7], myEdges[9], myEdges[11]);

		myEdges.clear();
		myTriangles.clear();
	}
	
	// Fill the simplices
	fillSimplices();
}

void
Cartosphere::TriangularMesh::reportAreas()
{
	FLP area = 0;
	FLP totalArea = 0;

	for (auto &triangle : listofSimplices)
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

void
Cartosphere::TriangularMesh::fillSimplices()
{
	listofSimplices.clear();

	unsigned pointIndex[6];
	for (unsigned index = 0; index < listofTriangles.size(); ++index)
	{
		auto& triplet = listofTriangles[index];
		// Extract the 3 vertices with duplication
		pointIndex[0] = listofEdges[std::get<0>(triplet).first].first;
		pointIndex[1] = listofEdges[std::get<0>(triplet).first].second;
		if (!std::get<0>(triplet).second)
			std::swap(pointIndex[0], pointIndex[1]);
		pointIndex[2] = listofEdges[std::get<1>(triplet).first].first;
		pointIndex[3] = listofEdges[std::get<1>(triplet).first].second;
		if (!std::get<1>(triplet).second)
			std::swap(pointIndex[2], pointIndex[3]);
		pointIndex[4] = listofEdges[std::get<2>(triplet).first].first;
		pointIndex[5] = listofEdges[std::get<2>(triplet).first].second;
		if (!std::get<2>(triplet).second)
			std::swap(pointIndex[4], pointIndex[5]);

		// Check if the edges form a simplex
		if (pointIndex[0] == pointIndex[1] || pointIndex[1] != pointIndex[2] ||
			pointIndex[2] == pointIndex[3] || pointIndex[3] != pointIndex[4] ||
			pointIndex[4] == pointIndex[5] || pointIndex[5] != pointIndex[0])
		{
			std::string message;
			{
				std::stringstream sst;
				sst << "Error in Face #" << index
					<< ": Edges do not form a valid simplex";
				message = sst.str();
			}
			listofMessages.push_back(message);
			flagofParsing = false;
			return;
		}

		// Push valid simplex with coordinates
		listofSimplices.emplace_back(
			listofPoints[pointIndex[0]],
			listofPoints[pointIndex[2]],
			listofPoints[pointIndex[4]]);
	}
	flagofParsing = true;
}
