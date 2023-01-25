
#ifndef __SHAPEFILE_HPP__
#define __SHAPEFILE_HPP__

#include <string>
using std::string;

namespace Cartosphere
{
    // ESRI-compliant ShapeFile Class
    class ShapeFile
    {
    public:
        // Values for shape type, p. 4
        enum ShapeType {
            NullShape = 0,
            Point = 1,
            PolyLine = 3,
            Polygon = 5,
            MultiPoint = 8,
            PointZ = 11,
            PolyLineZ = 13,
            PolygonZ = 15,
            MultiPointZ = 18,
            PointM = 21,
            PolyLineM = 23,
            PolygonM = 25,
            MultiPointM = 28,
            MultiPatch = 31,
        };
        
    public:
        // Opens shape file
        bool open(const std::string &folder, std::string &error);
        
        // Count number of shapes
        int count() const { return fileCode; }
        
    protected:
        // File specification
        int fileCode, fileLength, version;
        
        // Shape type (all non-null shapes are identical)
        ShapeType shapeType;
        
        // Bounding box (along x- and y-axes)
        int xMin, yMin, xMax, yMax;
        
        // Bounding box (along z- and m-axes)
        int zMin, zMax, mMin, mMax;
    };
}

#endif // !__SHAPEFILE_HPP__
