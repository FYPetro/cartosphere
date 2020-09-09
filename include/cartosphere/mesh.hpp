
#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <utility>
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
		// Default Constructor
		TriangularMesh() = default;
		// Construct triangular mesh from file
		TriangularMesh(std::string const &path) { load(path); }

	public:
		// Get file load state
		bool isLoaded() const { return flagofLoading; }
		// Get file load messages
		std::vector<std::string> getMessages() const { return listofMessages; }

	public:
		// Clear file
		void clear();
		// Load file from path
		bool load(std::string const &path);
		// Refine the mesh
		TriangularMesh refine(unsigned level = 1) const;
		// Report the area of each triangle
		void reportAreas();

	protected:
		// List of points
		std::vector<Point> listofPoints;
		// List of edges
		std::vector<std::pair<unsigned, unsigned>> listofEdges;
		// List of triangles
		std::vector<Triangle> listofTriangles;

	private:
		// File load flag
		bool flagofLoading = false;
		// File load message
		std::vector<std::string> listofMessages;
	};
}

#endif // !__MESH_HPP__
