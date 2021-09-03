
#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <utility>
#include <functional>
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
		Preimage(const Preimage& that) = default;
		// Copy assignment
		Preimage& operator=(const Preimage& that)
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
		Image(const FL3& f) : FL3(f) {}
		// Constructor from coordinates
		Image(FLP x, FLP y, FLP z) : FL3(x, y, z) {}
		// Copy constructor
		Image(const Image& that) = default;
		// Assignment constructor
		Image& operator=(const Image& that)
		{
			if (this != &that)
			{
				FL3::operator=(that);
			}
			return *this;
		}

	public:
		// Convert to vector
		FL3 to_vector() const { return static_cast<FL3>(*this); }
		// Obtain preimage
		Preimage to_preimage() const;
		// Obtain a unit vector
		Image to_unit_vector() const {
			return Image(static_cast<FL3>(*this).toUnitVector());
		}

	public:
		// Implicit conversion to FL3
		operator FL3() const { return static_cast<FL3>(*this); }
		// Calculate the spherical distance between two points
		friend FLP distance(const Image& a, const Image& b);
		// Calculate the angle spanned by three points
		friend FLP angle(const Image& a, const Image& b, const Image& c);
	};

	// A point keeps track of both coordinates
	class Point
	{
	public:
		// Default constructor, constructing an inconsistent point
		Point() : _preimage(), _image() {}
		// Copy constructor.
		Point(const Point& that) = default;
		// Construct from Preimage
		Point(const Preimage& preimage) :
			_preimage(preimage), _image(preimage.toImage()) {}
		// Construct from Image
		Point(const Image& image) :
			_preimage(image.to_preimage()), _image(image) {}

	public:
		// Get preimage
		inline Preimage const& preimage() const { return _preimage; }
		// Get image
		inline Image const& image() const { return _image; }
		// Obtain polar angle
		inline FLP p() const { return _preimage.p; }
		// Obtain azimuthal angle
		inline FLP a() const { return _preimage.a; }
		// Obtain x-coordinate
		inline FLP x() const { return _image.x; }
		// Obtain y-coordinate
		inline FLP y() const { return _image.y; }
		// Obtain z-coordinate
		inline FLP z() const { return _image.z; }

	public:
		// Flip to the antipode
		void flip() { _image.x *= -1; _image.y *= -1; _image.z *= -1; _populate_preimage(); }
		// Check if point is antipodal to another point
		bool isAntipodalTo(const Point& p) const;
		// Check if point is valid
		bool isValid() const;
		// Set preimage and fill image
		void set(const Preimage& preimage)
		{
			_preimage = preimage;
			_populate_image();
		}
		// Set image and fill preimage
		void set(const Image& image)
		{
			_image = image;
			_populate_preimage();
		}

	public:
		// Set preimage and fill image
		Point& operator=(const Preimage& that) { set(that); return *this; }
		// Set image and fill preimage
		Point& operator=(const Image& that) { set(that); return *this; }
		// Implicit conversion to preimage
		operator Preimage() const { return preimage(); }
		// Implicit conversion to image
		operator Image() const { return image(); }
		// Obtain the azimuth of a different given point relative to this point
		// Returns values in the interval (-pi,pi)
		FLP azimuth(const Point& other) const;
		// Obtain distance
		friend FLP distance(const Point& a, const Point& b)
		{
			return distance(a.image(), b.image());
		}
		// Calculate the angle spanned by three points
		friend FLP angle(const Point& a, const Point& b, const Point& c)
		{
			return angle(a.image(), b.image(), c.image());
		}
		// Obtain midpoint
		friend Point midpoint(const Point& a, const Point& b);

	private:
		// The preimage
		Preimage _preimage;
		// The image
		Image _image;

	protected:
		// Fill preimage from image
		void _populate_preimage() { _preimage = _image.to_preimage(); }
		// Fill image from preimage
		void _populate_image() { _image = _preimage.toImage(); }
	};

	typedef std::function<FLP(const Point&)> Function;

	// Representation of a directional minor arc and its local coordinate system
	class Arc
	{
	public:
		// Default Constructor
		Arc() {}
		// Constructor from two points
		Arc(const Point& A, const Point& B) : _a(A), _b(B) { _populate(); }

	public:
		// Return the angle spanned by the arc
		FLP span() const { return distance(_a.image(), _b.image()); }
		// Return the length of the arc
		FLP length() const { return distance(_a.image(), _b.image()); }
		// Return the pole of the arc
		Image pole() const { return Image(_n); }
		// Return the (t,0) local coordinates
		Image local(FLP u) const
		{
			return _a.image() * cos(u) + _c * sin(u);
		}
		// Return the (t,n) local coordinates
		Image local(FLP u, FLP v) const
		{
			Arc arc(local(u), Point(Image(_n)));
			return arc.local(v);
		}

	protected:
		// Populate the auxiliary vectors
		void _populate()
		{
			_n = normalize(cross(_a.image(), _b.image()));
			_c = cross(_n, _a.image());
		}

	protected:
		// Start and end points for the arc
		Point _a, _b;
		// The normal vector
		FL3 _n;
		// The point pi/2 radians away in the direction of the arc
		FL3 _c;
	};

	// A (spherical) triangle
	class Triangle
	{
	public:
		// A selection of quadrature rules for integration
		enum class Integrator {
			Centroid, ThreeVertices, Simpsons,
			Refinement1,
			Refinement2,
			Refinement3,
			Refinement4,
			Refinement5,
			Refinement6,
			Refinement7,
			Refinement8,
			Refinement9,
			Refinement10
		};
		static const Integrator DefaultIntegrator = Integrator::Refinement3;

	public:
		// Vertices
		Point A, B, C;

	public:
		// Default constructor
		Triangle() : A(), B(), C() {}
		// Construct with three points
		Triangle(const Point& P, const Point& Q, const Point& R) :
			A(P), B(Q), C(R) {}
		// Copy constructor
		Triangle(const Triangle& that) = default;

	public:
		// Obtain the orientation of the spherical triangle
		int orientation() const;
		// Calculate the area as a non-oriented spherical triangle
		// Returns an output in (0,pi)
		FLP area() const;
		// Returns the area of the oriented spherical triangle
		// Returns an output in (-pi,pi)
		FLP areaOriented() const { return orientation() * area(); }
		// Calculate the area as a Euclidean triangle
		// Returns a positive output
		FLP areaEuclidean() const;
		// Calculate the location of the center of mass
		Point centroid() const;
		// Obtain a finite element
		Function element(size_t index) const;
		// Numerically integrate a scalar function
		FLP integrate(const Function& f, Integrator intr) const;
	};

	// A (spherical) polygon
	class Polygon
	{
	public:
		// Default constructor
		Polygon() {}
		// Construct from a list of points
		Polygon(const std::vector<Point>& points) : _V(points) {}

	public:
		// Calculate the area as a spherical polygon
		FLP area() const;

	protected:
		// Input data: List of vertices
		std::vector<Point> _V;
	};

	// A triangular mesh
	class TriangularMesh
	{
	public:
		// Pair of indices to points
		typedef std::pair<size_t, size_t> UndirectedEdge;
		// Edge index and its orientation
		typedef std::pair<size_t, bool> DirectedEdge;
		// Triplet of indices with their orientations
		typedef std::tuple<DirectedEdge, DirectedEdge, DirectedEdge> DirectedEdgeTriplet;
		// Mesh statistics
		typedef struct
		{
			// The number of points
			size_t V;
			// The number of edges
			size_t E;
			// The number of triangles
			size_t F;
			// The area of the largest triangle
			FLP areaElementMax;
			// The area of the smallest triangle
			FLP areaElementMin;
			// The area ratio of the largest and smallest triangles
			FLP areaElementDisparity;
		} Stats;
		// Integrator rules
		enum class Quadrature {
			AreaWeighted
		};

	public:
		// Default Constructor
		TriangularMesh() = default;
		// Construct triangular mesh with a single triangle
		TriangularMesh(const Triangle& t) { load(t); }
		// Construct triangular mesh from file
		TriangularMesh(const std::string& path) { load(path); }

	public:
		// Get file load state
		bool isReady() const { return _bParseSuccess; }
		// Get file load messages
		std::vector<std::string> getMessages() const { return _vInfo; }
		// Get total area (spherical)
		FLP area() const;
		// Get total area (Euclidean)
		FLP areaEuclidean() const;

	public:
		// Clear file
		void clear();
		// Load triangle
		bool load(const Triangle& t);
		// Load file from path
		bool load(const std::string& path);
		// Save mesh to file
		bool save(const std::string& path) const;
		// Export mesh to OBJ format
		bool format(const std::string& path,
			const std::vector<FLP> &values = std::vector<FLP>()) const;
		// Export colored polyhedral object to OBJ format
		bool formatPoly(const std::string& path,
			const std::vector<FLP>& values) const;
		// Apply mid-point refinement to the mesh
		void refine();
		// Refine the mesh to a certain number of divisions
		void refine(size_t division);
		// Report the area of each triangle
		void reportAreas();
		// Numerically integrate a scalar function
		FLP integrate(const Function& f,
			Quadrature rule = Quadrature::AreaWeighted,
			Triangle::Integrator intr = Triangle::DefaultIntegrator) const;
		// Numerically integrate function values at vertices
		FLP integrate(const std::vector<FLP>& values) const;
		// Interpolate
		FLP interpolate(const std::vector<FLP>& values, const Point& point) const;
		// FEM: Generate inner products of the gradients of finite elements
		void fill(Matrix& A,
			Triangle::Integrator intr = Triangle::DefaultIntegrator) const;
		// FEM: Generate inner products of finite elements and their gradients
		void fill(Matrix& A, Matrix& M,
			Triangle::Integrator intr = Triangle::DefaultIntegrator) const;
		// FEM: Discretize an external force parametrized by x, y, and z.
		void fill(Vector& b, Function f,
			Triangle::Integrator intr = Triangle::DefaultIntegrator) const;
		// Generate statistics
		Stats statistics() const;
		// Return a list of vertices
		std::vector<Point> vertices() const;

	private:
		// Refresh redundant states
		void _populate();

	private:
		// Input data: List of points
		std::vector<Point> _V;
		// Input data: List of edges (using point indices)
		std::vector<UndirectedEdge> _E;
		// Input data: List of faces (each a triplet of directed edge indices)
		std::vector<DirectedEdgeTriplet> _F;
		// Redundant state: List of triangles (owning coordinates)
		std::vector<Triangle> _vt;
		// Redundant state: List of edges sharing a vertex
		std::vector<std::vector<size_t>> _VE;
		// Redundant state: List of faces sharing a vertex
		std::vector<std::vector<size_t>> _VF;
		// Redundant state: List of vertices in each face
		std::vector<std::vector<size_t>> _FV;

	private:
		// File load flag
		bool _bLoadSuccess = false;
		// File parse flag
		bool _bParseSuccess = false;
		// Messages
		std::vector<std::string> _vInfo;
	};
}

#endif // !__MESH_HPP__
