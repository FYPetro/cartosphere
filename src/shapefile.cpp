
#include "cartosphere/shapefile.hpp"
using namespace Cartosphere;

#include "cartosphere/utility.hpp"

#include <filesystem>
namespace fs = std::filesystem;
using fs::path;

#include <fstream>
using std::ifstream;

#include <iostream>
using std::cout;

bool
ShapeFile::open(const std::string &folder, std::string &error)
{
    // Check if the folder exists.
    path directory{folder};
    if (!fs::exists(directory))
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
            // cout << "New Polygon: " << shapes.size() << "\n";

            ptr = std::make_shared<Polygon>();
            Polygon& polygon = *static_cast<Polygon*>(ptr.get());
            polygon.type = recordShapeType;

            // cout << "Type: " << polygon.type << "\n";

            file.read(reinterpret_cast<char*>(floats), 4 * sizeof(double));
            polygon.box[0] = floats[0];
            polygon.box[1] = floats[1];
            polygon.box[2] = floats[2];
            polygon.box[3] = floats[3];

            // cout << "Box:";
            // for (int i = 0; i < 4; ++i)
            // {
            //     cout << " " << polygon.box[i];
            // }
            // cout << "\n";

            file.read(reinterpret_cast<char*>(integers), 2 * sizeof(int));
            polygon.numParts = integers[0];
            polygon.numPoints = integers[1];

            // cout << "Parts: " << polygon.numParts << "\n";
            // cout << "Points: " << polygon.numPoints << "\n";

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

        // cout << parsedLength << "/" << fileLength <<  "\n";
    }
    
    file.close();
    
    // ================ //
    // Parse Index File //
    file.open(directory/idxFile, std::ios::binary);
    if (!file.is_open())
    {
        error = "The index file " + idxFile + " does not exist.";
        return false;
    }
    
    file.open(directory/dbFile, std::ios::binary);
    if (!file.is_open())
    {
        error = "The dBASE file " + dbFile + " does not exist.";
        return false;
    }
    
    return true;
}
