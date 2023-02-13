
#ifndef __SHAPEFILE_HPP__
#define __SHAPEFILE_HPP__

#include "cartosphere/mesh.hpp"

namespace Cartosphere
{
	// ESRI-compliant ShapeFile Class
	class ShapeFile
	{
	public:
		// Values for shape type, p. 4
		enum ShapeType {
			NullShapeType = 0,
			PointType = 1,
			PolyLineType = 3,
			PolygonType = 5,
			MultiPointType = 8,
			PointZType = 11,
			PolyLineZType = 13,
			PolygonZType = 15,
			MultiPointZType = 18,
			PointMType = 21,
			PolyLineMType = 23,
			PolygonMType = 25,
			MultiPointMType = 28,
			MultiPatchType = 31,
		};

		// Null Shape (0)
		class Shape
		{
		public:
			// Shape type
			ShapeType type;

		public:
			// Print to matlab
			virtual string to_matlab() const = 0;
		};

		// Alias for shared pointers to Shape
		typedef shared_ptr<Shape> MyShapePtr;

		// Type Point (1)
		class Point : public Shape
		{
		public:
			// X coordinate
			double x;

			// Y coordinate
			double y;

		public:
			// Print to matlab
			string to_matlab() const;
		};

		// Type Polygon (5)
		class Polygon : public Shape
		{
		public:
			// Bounding box
			double box[4];

			// Number of parts
			int numParts;

			// Total number of points
			int numPoints;

			// Index to first point in part
			vector<int> parts;

			// Points for all parts
			vector<Point> points;

		public:
			// Output to matlab
			virtual string to_matlab() const;
		};

	public:
		// Opens shape file
		bool open(const string &folder, string &error);

		// Outputs shape file as a matlab file
		void to_matlab(const string &path) const;

		// Offloads all points into a single vector
		vector<Cartosphere::Point> gather() const;

		// Count number of shapes
		size_t count() const { return shapes.size(); }

	protected:
		// File specification
		int fileCode, fileLength, version;

		// Shape type (all non-null shapes are identical)
		ShapeType shapeType;

		// Bounding box (along x- and y-axes)
		double xMin, yMin, xMax, yMax;

		// Bounding box (along z- and m-axes)
		double zMin, zMax, mMin, mMax;

		// Shapes parsed
		vector<MyShapePtr> shapes;
	};
}

#endif // !__SHAPEFILE_HPP__
