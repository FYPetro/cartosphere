
#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <vector>
#include <utility>
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
		Point &operator=(Preimage const &that) { set(that); return *this; }
		Point &operator=(Image const &that) { set(that); return *this; }
		// Obtain polar angle
		// FLP p() const { return pi.p; }
		// Obtain azimuthal angle
		// FLP a() const { return im.a; }
		// Obtain x-coordinate
		// FLP x() const { return im.x; }
		// Obtain y-coordinate
		// FLP y() const { return im.y; }
		// Obtain z-coordinate
		// FLP z() const { return im.z; }

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
		TriangularMesh(std::string const &path) {}

	protected:
		// List of points
		std::vector<Point> listofPoints;
		// Point membership of triangles
		std::vector<UI3> listofMembers;
		// List of triangles
		std::vector<Triangle> listofTriangles;

		TriangularMesh refine(unsigned level = 1) const;
	};
}

#endif // !__MESH_HPP__
