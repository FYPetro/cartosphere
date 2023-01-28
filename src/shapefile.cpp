
#include "cartosphere/shapefile.hpp"
using Cartosphere::ShapeFile;
using std::shared_ptr;
using std::string;
using std::vector;

#include "cartosphere/utility.hpp"

#include <filesystem>
namespace filesystem = std::filesystem;

string
ShapeFile::Point::to_matlab() const
{
	stringstream sst;

	return sst.str();
}

string
ShapeFile::Polygon::to_matlab() const
{
	stringstream sst;

	vector<int> delimiters = parts;
	delimiters.push_back((int)points.size());

	// for (int k = 0; k <= numParts; ++k)
	// {
	// 	std::cout << delimiters[k] << ",";
	// }
	// std::cout << "\n";

	// Create a column cell array to store the polylines
	sst << "temp = cell(" << numParts << ",1);\n";
	for (int k = 0; k < numParts; ++k)
	{
		sst << "temp{" << (k + 1) << "} = [...\n";
		for (int i = delimiters[k]; i < delimiters[k + 1]; ++i)
		{
			const auto& point = points[i];
			sst << "	" << point.x << " " << point.y << ";\n";
		}
		sst << "];\n";
	}

	return sst.str();
}

bool
ShapeFile::open(const string &folder, string &error)
{
	// Check if the folder exists.
	filesystem::path directory{folder};
	if (!filesystem::exists(directory))
	{
		error = "The folder does not exist.";
		return false;
	}

	// Prepare to read shape, index, and dBASE files.
	ifstream file;
	string mainFile = folder + ".shp";
	string idxFile = folder + ".shx";
	string dbFile = folder + ".dbf";

	// =============== //
	// Parse Main File //
	file.open(directory/mainFile, std::ios::binary);
	if (!file.is_open())
	{
		error = "The main file " + mainFile + " does not exist.";
		return false;
	}

	// Parse Main File Header
	int integers[9]{};
	for (int i = 0; i < 9 && file.good(); ++i)
	{
		file.read(reinterpret_cast<char *>(integers + i), sizeof(int));
		// First seven bytes are big-endian.
		if (i < 7)
		{
			endswap(integers + i);
		}
	}
	if (!file.good())
	{
		error = "Main file header error: integers.";
		return false;
	}
	fileCode = integers[0];
	fileLength = integers[6];
	version = integers[7];
	shapeType = (ShapeType)integers[8];

	double floats[8]{};
	for (int i = 0; i < 8 && file.good(); ++i)
	{
		floats[i] = 0;
		file.read(reinterpret_cast<char *>(floats + i), sizeof(double));
	}
	if (!file.good())
	{
		error = "Main file header error: floats.";
		return false;
	}
	xMin = floats[0]; yMin = floats[1];
	xMax = floats[2]; yMax = floats[3];
	zMin = floats[4]; zMax = floats[5];
	mMin = floats[6]; mMax = floats[7];

	// Parse Records (2 bytes counts as 1 word)
	int parsedLength = 100 / 2;
	while (parsedLength < fileLength)
	{
		// Read Record Header
		file.read(reinterpret_cast<char*>(integers), 2*sizeof(int));
		endswap(integers);
		endswap(integers + 1);
		int recordNumber = integers[0];
		int contentLength = integers[1];

		// Read Record Contents
		file.read(reinterpret_cast<char*>(integers), sizeof(int));
		ShapeType recordShapeType = (ShapeType)integers[0];

		MyShapePtr ptr;
		switch (shapeType)
		{
		case NullShapeType:
			break;

		case PolygonType:
		{
			ptr = std::make_shared<Polygon>();
			Polygon& polygon = *static_cast<Polygon*>(ptr.get());
			polygon.type = recordShapeType;

			file.read(reinterpret_cast<char*>(floats), 4 * sizeof(double));
			polygon.box[0] = floats[0];
			polygon.box[1] = floats[1];
			polygon.box[2] = floats[2];
			polygon.box[3] = floats[3];

			file.read(reinterpret_cast<char*>(integers), 2 * sizeof(int));
			polygon.numParts = integers[0];
			polygon.numPoints = integers[1];

			for (int i = 0; i < polygon.numParts; ++i)
			{
				file.read(reinterpret_cast<char*>(integers), sizeof(int));
				polygon.parts.push_back(integers[0]);
			}

			Point point;
			for (int i = 0; i < polygon.numPoints; ++i)
			{
				file.read(reinterpret_cast<char*>(floats), 2*sizeof(double));
				point.x = floats[0];
				point.y = floats[1];
				polygon.points.push_back(point);
			}
		}
		break;

		default:
		{
			error = "Not implemented.";
			return false;
		}
		}

		shapes.push_back(ptr);
		parsedLength += 4 + contentLength;

		// std::cout << parsedLength << "/" << fileLength <<  "\n";
	}
	file.close();

	// Parse Index File //
	file.open(directory/idxFile, std::ios::binary);
	if (!file.is_open())
	{
		error = "The index file " + idxFile + " does not exist.";
		return false;
	}
	file.close();

	file.open(directory/dbFile, std::ios::binary);
	if (!file.is_open())
	{
		error = "The dBASE file " + dbFile + " does not exist.";
		return false;
	}
	file.close();

	return true;
}

void
ShapeFile::to_matlab(const string &outPath) const
{
	ofstream ofs;

	// Print points of all shapes in *.m_data.m
	auto outStem = filesystem::path{outPath}.stem();
	auto outDataPath = outStem.filename().string() + "_data.m";
	ofs.open(outDataPath);
	ofs << "%% Shapes\n";
	for (size_t i = 0; i < shapes.size(); ++i)
	{
		ofs << "% Shape " << (i + 1) << "\n"
			<< shapes[i]->to_matlab()
			<< "shape{" << (i + 1) << "} = temp;\n";
	}
	ofs << "clear temp;\n";
	ofs.close();

	// Plot all shapes in *.m
	ofs.open(outPath);
	ofs << "%% Load data\n"
		<< outStem.string() << "_data\n\n"
		<< "%% Plot shape\n"
		<< "for i = 1:length(shape)\n"
		<< "\t" << "for k = 1:length(shape{i})\n"
		<< "\t\t" << "plot(shape{i}{k}(:,1),shape{i}{k}(:,2));\n"
		<< "\t\t" << "hold on\n"
		<< "\t"	<< "end\n"
		<< "end\n"
		<< "hold off\n"
		<< "axis equal\n\n";

	// Provide space for user code
	ofs << "%% Write your code here\n";
	ofs.close();
}
