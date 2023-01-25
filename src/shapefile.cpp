
#include "cartosphere/shapefile.hpp"
using namespace Cartosphere;

#include "cartosphere/utility.hpp"

#include <filesystem>
namespace filesystem = std::filesystem;
using filesystem::path;

#include <fstream>
using std::ifstream;

bool
ShapeFile::open(const std::string &folder, std::string &error)
{
    // Check if the folder exists.
    path directory{folder};
    if (!filesystem::exists(directory))
    {
        error = "The folder does not exist.";
        return false;
    }
    
    // Prepare to read shape, index, and dBASE files.
    ifstream file;
    string mainPath = folder + ".shp";
    string idxPath = folder + ".shx";
    string dbPath = folder + ".dbf";
    
    // =============== //
    // Parse Main File //
    file.open(directory/mainPath, std::ios::binary);
    if (!file.is_open())
    {
        error = "The main file " + mainPath + " does not exist.";
        return false;
    }
    
    // Parse Main File Header
    int ints[9];
    for (int i = 0; i < 9 && file.good(); ++i)
    {
        file.read(reinterpret_cast<char *>(ints + i), sizeof(int));
        // First seven bytes are big-endian.
        if (i < 7)
        {
            endswap(ints + i);
        }
    }
    if (!file.good())
    {
        error = "Header error.";
        return false;
    }
    fileCode = ints[0];
    fileLength = ints[6];
    version = ints[7];
    shapeType = (ShapeType)ints[8];
    
    double floats[8];
    for (int i = 0; i < 8 && file.good(); ++i)
    {
        file.read(reinterpret_cast<char *>(floats + i), sizeof(double));
    }
    if (!file.good())
    {
        error = "Header error.";
        return false;
    }
    xMin = floats[0]; yMin = floats[1];
    xMax = floats[2]; yMax = floats[3];
    zMin = floats[4]; zMax = floats[5];
    mMin = floats[6]; mMax = floats[7];
    
    // Parse Records
    int parsedLength = 100;
    
    file.close();
    
    // ================ //
    // Parse Index File //
    file.open(directory/idxPath, std::ios::binary);
    if (!file.is_open())
    {
        error = "The index file " + idxPath + " does not exist.";
        return false;
    }
    
    file.open(directory/dbPath, std::ios::binary);
    if (!file.is_open())
    {
        error = "The dBASE file " + dbPath + " does not exist.";
        return false;
    }
    
    return true;
}
