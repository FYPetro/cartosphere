
#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <utility>
#include <tuple>
#include <vector>
#include <string>

#include "cartosphere/utility.hpp"

namespace Cartosphere
{
	class Image;

	// Coordinates of a point using spherical parametrization
	class Preimage
	{
	public:
		// Polar angle
		FLP p;
		// Azimuthal
		FLP a;

	public:
		// Default constructor
		Preimage() : p(0), a(0) {}
		// Constructor from two angles
		Preimage(FLP polar, FLP azimuthal) : p(polar), a(azimuthal) {}
		// Copy constructor
		Preimage(Preimage const &that) = default;
		// Copy assignment
		Preimage &operator=(Preimage const &that)
		{
			if (this != &that)
			{
				p = that.p;
				a = that.a;
			}
			return *this;
		}

	public:
		// Obtain image
		Image toImage() const;
	};

	// Coordinates of a point embedded in the 3-space
	class Image : public FL3
	{
	public:
		// Default Constructor
		Image() {}
		// Construct from FL3
		Image(FL3 const &f) : FL3(f) {}
		// Constructor from coordinates
		Image(FLP x, FLP y, FLP z) : FL3(x, y, z) {}
		// Copy constructor
		Image(Image const &that) = default;

	public:
		// Calculate the spherical distance between two points
		friend FLP distance(Image const &a, Image const &b)
		{
			Preimage pa = a.toPreimage();
			Preimage pb = b.toPreimage();
			FLP value = cos(pa.p) * cos(pb.p) +
				sin(pa.p) * sin(pb.p) * cos(pa.a - pb.a);
			return acos(value);
		}

	public:
		// Obtain preimage
		Preimage toPreimage() const;
		// Assignment constructor
		Image &operator=(Image const &that)
		{
			if (this != &that)
				FL3::operator=(that);
			return *this;
		}
	};

	// A point keeps track of both coordinates
	class Point
	{
	public:
		// Default constructor, constructing an inconsistent point
		Point() : pi(), im() {}
		// Copy constructor.
		Point(Point const &that) = default;
		// Construct from Preimage
		Point(Preimage const &preimage) : pi(preimage), im(preimage.toImage()) {}
		// Construct from Image
		Point(Image const &image) : pi(image.toPreimage()), im(image) {}

	public:
		// Get preimage
		Preimage const &preimage() const { return pi; }
		// Get image
		Image const &image() const { return im; }
		// Obtain polar angle
		FLP p() const { return pi.p; }
		// Obtain azimuthal angle
		FLP a() const { return pi.a; }
		// Obtain x-coordinate
		FLP x() const { return im.x; }
		// Obtain y-coordinate
		FLP y() const { return im.y; }
		// Obtain z-coordinate
		FLP z() const { return im.z; }

	public:
		// Check if point is antipodal to another point
		bool isAntipodalTo(Point const &p) const
		{
			return this != &p
				&& this->x() + p.x() == 0
				&& this->y() + p.y() == 0
				&& this->z() + p.z() == 0;
		}
		// Check if point is valid
		bool isValid() const
		{
			return p() == 0 && a() == 0 && x() == 0 && y() == 0 && z() == 0;
		}
		// Set preimage and fill image
		void set(Preimage const &preimage)
		{
			pi = preimage;
			populateImage();
		}
		// Set image and fill preimage
		void set(Image const &image)
		{
			im = image;
			populatePreimage();
		}

	public:
		// Set preimage and fill image
		Point &operator=(Preimage const &that) { set(that); return *this; }
		// Set image and fill preimage
		Point &operator=(Image const &that) { set(that); return *this; }
		// Obtain midpoint
		friend Point midpoint(Point const &a, Point const &b);

	private:
		// The preimage
		Preimage pi;
		// The image
		Image im;

	protected:
		// Fill preimage from image
		void populatePreimage()	{ pi = im.toPreimage(); }
		// Fill image from preimage
		void populateImage() { im = pi.toImage(); }
	};

	// Representation of a directional minor arc and its local coordinate system
	class Arc
	{
	public:
		// Default Constructor
		Arc() {}
		// Constructor from two points
		Arc(Point const &A, Point const &B) : a(A), b(B) { fill(); }

	public:
		// Return the angle spanned by the arc
		FLP span() const { return distance(a.image(), b.image()); }
		// Return the length of the arc
		FLP length() const { return distance(a.image(), b.image()); }
		// Return the (t,0) local coordinates
		Point local(FLP u) const
		{
			return Point(a.image() * cos(u) + c * sin(u));
		}
		// Return the (t,n) local coordinates
		Point local(FLP u, FLP v) const
		{
			Arc arc(local(u), Point(Image(n)));
			return arc.local(v);
		}

	protected:
		// Populate the auxiliary vectors
		void fill()
		{
			n = normalize(cross(a.image(), b.image()));
			c = cross(n, a.image());
		}

	protected:
		// Start and end points for the arc
		Point a, b;
		// The normal vector
		FL3 n;
		// The point pi/2 radians away in the direction of the arc
		FL3 c;
	};

	// A (spherical) triangle
	class Triangle
	{
	public:
		// Vertices
		Point A, B, C;

	public:
		// Default constructor
		Triangle() : A(), B(), C() {}
		// Construct with three points
		Triangle(Point const &P, Point const &Q, Point const &R) : A(P), B(Q), C(R) {}
		// Copy constructor
		Triangle(Triangle const &that) = default;

	public:
		// Calculate the area as a spherical triangle
		FLP area() const;
		// Calculate the area as a Euclidean triangle
		FLP areaEuclidean() const;
	};

	// A triangular mesh
	class TriangularMesh
	{
	public:
		// Pair of indices to points
		typedef std::pair<unsigned, unsigned> PointIndexPair;
		// Edge index and its orientation
		typedef std::pair<unsigned, bool> EdgeIndex;
		// Triplet of indices with their orientations
		typedef std::tuple<EdgeIndex, EdgeIndex, EdgeIndex> EdgeIndexTriplet;
		// Memory statistics
		typedef struct
		{
			unsigned nPoint;
			unsigned nEdge;
			unsigned nTriangle;
			FLP areaElementMax;
			FLP areaElementMin;
			FLP areaElementDisparity;
		} Stats;

	public:
		// Default Constructor
		TriangularMesh() = default;
		// Construct triangular mesh from file
		TriangularMesh(std::string const &path) { load(path); }

	public:
		// Get file load state
		bool isReady() const { return flagofParsing; }
		// Get file load messages
		std::vector<std::string> getMessages() const { return listofMessages; }
		// Get total area (spherical)
		FLP area() const;
		// Get total area (Euclidean)
		FLP areaEuclidean() const;

	public:
		// Clear file
		void clear();
		// Load file from path
		bool load(std::string const &path);
		// Save mesh to file
		bool save(std::string const &path) const;
		// Export mesh to OBJ format
		bool format(std::string const &path) const;
		// Refine the mesh
		void refine();
		// Report the area of each triangle
		void reportAreas();
		// Generate statistics
		Stats statistics() const;

	private:
		// Fill simplicial complex
		void fillSimplices();

	private:
		// List of points
		std::vector<Point> listofPoints;
		// List of edges (using point indices)
		std::vector<PointIndexPair> listofEdges;
		// List of triangles (using edge indices)
		std::vector<EdgeIndexTriplet> listofTriangles;
		// List of triangles (owning coordinates)
		std::vector<Triangle> listofSimplices;

	private:
		// File load flag
		bool flagofLoading = false;
		// File parse flag
		bool flagofParsing = false;
		// Messages
		std::vector<std::string> listofMessages;
	};
}

#endif // !__MESH_HPP__
