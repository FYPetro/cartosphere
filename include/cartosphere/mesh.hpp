
#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <utility>
#include <functional>
#include <iterator>
#include <tuple>
#include <set>
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
		// Construct from Preimage parameters
		Point(FLP polar, FLP azimuthal) :
			Point(Preimage(polar, azimuthal)) {}
		// Construct from Image
		Point(const Image& image) :
			_preimage(image.to_preimage()), _image(image) {}
		// Construct from Image parameters
		Point(FLP x, FLP y, FLP z) :
			Point(Image(x, y, z)) {}

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
		// Move point along a tangent displacement
		void move(FL3 displacement);

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
		Arc(const Point& A, const Point& B) :
			_a(A), _b(B)
		{ _populate(); }
		// Constructor from a point and a direction vector
		Arc(const Point& A, FL3 u) :
			_a(A), _b(normalize(A.image() + u))
		{ _populate(); }

	public:
		// Return the signed distance of a point to the arc
		FLP distance(const Point& p) const
		{
			return M_PI_2 - acos(dot(_n, p));
		}
		// Return if a point is in the left hemisphere
		bool encloses(const Point& p) const
		{
			return dot(_n, p) >= 0;
		}
		// Return the angle spanned by the arc
		FLP span() const
		{
			return length();
		}
		// Returns a tangent vector at a local coordinate
		FL3 tangent(FLP u) const
		{
			return _a.image() * (-sin(u)) + _c * cos(u);
		}
		// Return the length of the arc (on S2)
		FLP chord() const
		{
			return (_a.image() - _b.image()).norm2();
		}
		// Return the length of the arc (on S2)
		FLP length() const
		{
			return atan2(_sin, _cos);
		}
		// Return the midpoint
		Point midpoint() const { return Point((_a.image() + _b.image()).normalize()); }
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
		// Rotate by arc
		FL3 rotate(FL3 vector) const
		{
			// Rodrigues' rotation formula.
			FL3 rotated = vector * _cos
				+ cross(_n, vector) * _sin
				+ _n * (dot(_n, vector) * (1 - _cos));
			return rotated;
		}

	protected:
		// Populate auxiliary states
		void _populate()
		{
			// Use axb to determine:
			// 1. a normal vector in the correct orientation
			// 2. the sine of the subtended angle
			// Note that a and b must be given as unit vectors.
			FL3 axb = cross(_a.image(), _b.image());
			_n = normalize(axb);
			_c = cross(_n, _a.image());
			_cos = dot(_a.image(), _b.image());
			_sin = axb.norm2();
		}

	protected:
		// Start and end points for the arc
		Point _a, _b;
		// The normal vector
		FL3 _n;
		// The point pi/2 radians away in the direction of the arc
		FL3 _c;
		// Trigs of angle
		FLP _cos, _sin;
	};

	// Representation of a spherical cap by its cente and radius
	class Cap
	{
	public:
		// Default constructor
		Cap() = default;
		// Constructor from apex and spherical radius
		Cap(const Point& apex, FLP radius) : _a(apex), _r(radius) {}

	public:
		// Circumscription calculation
		template <class ForIt>
		static Cap circumscribe(ForIt first, ForIt last)
		{
			// TODO
			return Cap();
		}

	public:
		// Check for point containment
		bool contains(const Point& p) const { return distance(_a, p) <= _r; }
		// Returns the apex
		Point apex() const { return _a; }
		// Returns the radius
		FLP radius() const { return _r; }
		// Combines with another cap
		//Cap combine(const Cap& other) const;

	protected:
		// Apex
		Point _a;
		// Radius
		FLP _r = 0;
	};

	// A (spherical) triangle
	class Triangle
	{
	public:
		// A selection of quadrature rules for integration
		enum class Integrator {
			Centroid,
			ThreeVertices,
			Simpsons,
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
		// Compute the barycentric coordinates
		FL3 barycentric(const Point& p) const;
		// Calculate the location of the center of mass
		Point centroid() const;
		// Check for point containment
		bool contains(const Point& p) const;
		// Calculate the diameter of the circumcircle
		FLP diameter() const;
		// Compute a circumcircle
		Cap circumcircle() const;
		// Obtain a finite element
		Function element(size_t index) const;
		// Obtain a gradient vector element
		FL3 gradient(size_t index) const;
		// Numerically integrate a scalar function
		FLP integrate(const Function& f, Integrator intr) const;
	};

	// A (spherical) polygon
	class Polygon
	{
	public:
		enum class Relation
		{
			Interior, Boundary, Outside
		};

	public:
		// Default constructor
		Polygon() {}

		// Construct from a list of points
		Polygon(const std::vector<Point>& points) : _V(points) {}

	public:
		// Calculate the area as a spherical polygon
		FLP area() const;

		// Check if a point is in the interior.
		Relation interior(const Point& point) const;

	public:
		// Emplace a point using spherical coordinates
		inline void emplace_back(FLP polar, FLP azimuth)
		{
			_V.emplace_back(polar, azimuth);
		}

		// Emplace a point using Cartesian coordinates
		inline void emplace_back(FLP x, FLP y, FLP z)
		{
			_V.emplace_back(x, y, z);
		}

		// Remove all points
		inline void clear()
		{
			_V.clear();
		}

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
			// The diameter of the largest circumcircle
			FLP diameterElementMax;
		} Stats;
		// Integrator rules
		enum class Quadrature {
			AreaWeighted,
			DualAreaWeighted, // Weight values by area of dual polygons
		};
		// SSTree
		class Tree
		{
		protected:
			// SSElem
			class Node
			{
			public:
				// Type: Node Pointer
				typedef std::shared_ptr<Node> Pointer;
				// Type: Vector of pointers type
				typedef std::vector<Pointer> Pointers;
				// Why not?
				static const size_t npos = (size_t)-1;

			protected:
				// Stores the parent
				Pointer _parent;
				// Stores the children
				Pointers _children;
				// Index of the triangle
				size_t _index = npos;
				// Number of nodes in the subtree (total number of children)
				size_t _subtree = 0;
				// The height above this node
				size_t _height = 0;
				// Update count without refresh values
				size_t _updates = 0;
				// The sum of squared distance
				FLP _variance = 0;
				// Circumcircle, containing the centroid and a radius (S2-metric)
				Cap _cap;

			public:
				// Make a node
				Node() = default;
				// Make a node from an index and a circumcircle
				Node(size_t index, const Cap& cap) : _index(index), _cap(cap) {}

			public:
				// Number of children
				size_t degree() const { return _children.size(); }
				// Height of node
				size_t height() const { return _height; }
				// Return parent
				Pointer parent() const { return _parent; }
				// Is it a leaf?
				bool isLeaf() const { return _children.empty(); }
				// Is it a root?
				bool isRoot() const { return _parent == nullptr; }
				// Obtain a cap
				const Cap& cap() const { return _cap; }
				// Obtain all children
				Pointers children() const { return _children; }
				// Obtain pointer to first children
				Pointers::const_iterator begin() const { return _children.begin(); }
				// Obtain pointer to first children
				Pointers::const_iterator end() const { return _children.end(); }

			public:
				// Update the children and the circumcircle
				void update(const Node::Pointers& nodes)
				{
					_children = nodes;
					std::vector<Cap> caps;
					std::transform(_children.begin(), _children.end(), std::back_inserter(caps),
						[](const Pointer& p) {
							return p->cap();
						}
					);
					_cap = Cap::circumscribe(caps.begin(), caps.end());
				}

			public:
				// Accommodate an entry
				void store(Node::Pointer n) { _children.push_back(n); }
			};
			// Minimum number of children
			static const size_t m = 2;
			// Maximum number of children
			static const size_t M = 80;
			// Number of points to omit when forcing a reinsertion
			static const size_t p = 32;

		public:
			// Is empty?
			bool empty() const;
			// Index a triangle based on a point
			size_t find(const Point& p) const;
			// Height of the tree
			size_t height() const { return _root.height(); }

		public:
			// Clear the tree
			void clear();
			// Build from triangles
			void build(const std::vector<Triangle>& vt)
			{
				// Clear the tree
				clear();
				// Insert the triangles into the tree one by one.
				for (size_t i = 0; i < vt.size(); ++i)
				{
					insert(i, vt[i].circumcircle());
				}
			}
			// Insert a circumcircle with its index
			void insert(size_t index, const Cap& cap)
			{
				_insert(std::make_shared<Tree::Node>(index, cap));
			}

		protected:
			// Algorithm Insert
			void _insert(Node::Pointer n, size_t level = 0);
			// Algorithm ChooseSubtree
			Node::Pointer _choose(Node::Pointer n, size_t level = 0);
			// Algorithm OverflowTreatment
			bool _overflow(Node::Pointer n, size_t level = 0);
			// Algorithm ReInsert
			void _reinsert(Node::Pointer n);
			// Algorithm Split
			void _split(Node::Pointer n);

		protected:
			// Root node
			Node _root;
			// Overflow calls
			std::vector<bool> _overflown;
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
		// Return barycentric coordinates and provide the triangle through reference
		FL3 barycentric(const Triangle& t) const;

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
		// Compute the gradient field at the vertices
		void set(const std::vector<FLP>& values);
		// Interpolate the function value of a point
		FLP interpolate(const Point& p) const;
		// Interpolate the gradient value at a point
		FL3 gradient(const Point& p) const;
		// Numerically integrate a scalar function
		FLP integrate(const Function& f,
			Quadrature rule = Quadrature::AreaWeighted,
			Triangle::Integrator intr = Triangle::DefaultIntegrator) const;
		// Return the index vertices of the triangles
		UI3 indexTriangleVertices(size_t triangleIndex) const
		{
			std::vector<size_t> V = _FV[triangleIndex];
			UI3 index;
			index.a = V[0];
			index.b = V[1];
			index.c = V[2];
			return index;
		}
		// Numerically integrate function values at vertices
		FLP integrate(const std::vector<FLP>& values,
			Quadrature rule = Quadrature::DualAreaWeighted,
			Triangle::Integrator intr = Triangle::DefaultIntegrator) const;
		// L2 error
		FLP lebesgue(const std::vector<FLP>& weights, const Function& func,
			Triangle::Integrator intr = Triangle::DefaultIntegrator) const;
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
		std::vector<Point> vertices() const
		{
			return _V;
		}

	private:
		// Refresh redundant states
		void _populate();
		// Compute the gradient given nodal values
		void _gradient(const std::vector<FLP>& a);
		// Lookup triangle index from a point
		size_t _lookup(const Point& p) const;

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
		// State: Nodal values
		std::vector<FLP> _a;
		// State: gradient vectors
		std::vector<FL3> _grad;
		// State: SS-tree for lookup
		// Tree _tree;

	private:
		// File load flag
		bool _bLoadSuccess = false;
		// File parse flag
		bool _bParseSuccess = false;
		// Messages
		std::vector<std::string> _vInfo;
	};

	// Parallel transport of a vector
	FL3 transport(const Point& from, const Point& to, const FL3& vector);
}

#endif // !__MESH_HPP__
