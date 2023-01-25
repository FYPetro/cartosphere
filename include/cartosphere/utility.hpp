
#ifndef __UTILITY_HPP__
#define __UTILITY_HPP__

using FLP = double;

#include <cmath>

#pragma warning(push)
#pragma warning(disable: 4819)
#include <Eigen/Core>

#include <Eigen/Sparse>

#include <Eigen/SparseCore>

#include <Eigen/IterativeLinearSolvers>
#pragma warning(pop)

using Matrix = Eigen::SparseMatrix<FLP, Eigen::RowMajor>;
using Entry = Eigen::Triplet<FLP>;
using Vector = Eigen::Matrix<FLP, Eigen::Dynamic, 1>;
using Solver = Eigen::BiCGSTAB<Matrix,Eigen::IncompleteLUT<Matrix::Scalar>>;

template<typename T>
constexpr FLP deg2rad(T angle) { return (angle * (FLP)M_PI / (FLP)180); }

template<typename T>
constexpr FLP rad2deg(T angle) { return (angle * (FLP)M_1_PI * (FLP)180); }

static const FLP EPS = std::numeric_limits<FLP>::epsilon();

class FL3
{
public:
	// Components
	FLP x, y, z;

public:
	// Default Constructor
	FL3() : x(0), y(0), z(0) {}
	// Construct from components
	FL3(FLP a, FLP b, FLP c) : x(a), y(b), z(c) {}
	// Copy constructor
	FL3(const FL3& that) = default;
	// Assignment operator
	FL3& operator=(const FL3& that)
	{
		if (this != &that)
		{
			x = that.x;
			y = that.y;
			z = that.z;
		}
		return *this;
	}

public:
	// Scalar multiplication
	FL3 operator*(FLP c) const { return FL3(x * c, y * c, z * c); }
	// Negation
	FL3 operator-() const { return FL3(-x, -y, -z); }
	// Add a vector
	FL3& operator+=(FL3 const& b)
	{
		this->x -= b.x;
		this->y -= b.y;
		this->z -= b.z;
		return *this;
	}
	// Scalar multiplication by constant
	FL3& operator*=(FLP b)
	{
		x *= b;
		y *= b;
		z *= b;
		return *this;
	}
	// Component-wise division by constant
	FL3& operator/=(FLP b)
	{
		x /= b;
		y /= b;
		z /= b;
		return *this;
	}
	// Component-wise addition
	friend FL3 operator+(const FL3& a, const FL3& b)
	{
		return FL3(a.x + b.x, a.y + b.y, a.z + b.z);
	}
	// Component-wise subtraction
	friend FL3 operator-(const FL3& a, const FL3& b)
	{
		return FL3(a.x - b.x, a.y - b.y, a.z - b.z);
	}
	// Scalar multiplication
	friend FL3 operator*(FLP c, const FL3& b)
	{
		return b * c;
	}
	// Component-wise division by constant
	friend FL3 operator/(const FL3& a, FLP b)
	{
		return FL3(a.x / b, a.y / b, a.z / b);
	}
	// Dot product
	friend FLP dot(const FL3& a, const FL3& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}
	// Cross product
	friend FL3 cross(const FL3& a, const FL3& b)
	{
		FLP x = a.y * b.z - a.z * b.y;
		FLP y = a.z * b.x - a.x * b.z;
		FLP z = a.x * b.y - a.y * b.x;
		return FL3(x, y, z);
	}
	// Triple product
	friend FLP triple(const FL3& a, const FL3& b, const FL3& c)
	{
		return a.x * b.y * c.z + a.y * b.z * c.x + a.z * b.x * c.y
			- (a.z * b.y * c.x + a.x * b.z * c.y + a.y * b.x * c.z);
	}
	// Normalize
	friend FL3 normalize(const FL3& a)
	{
		FL3 unit = a;
		return unit.normalize();
	}

public:
	// Returns the square of the 2-norm
	FLP norm2sq() const { return pow(x, 2) + pow(y, 2) + pow(z, 2); }
	// Returns the 2-norm
	FLP norm2() const { return sqrt(norm2sq()); }
	// Normalizes the vector
	FL3& normalize() { (*this) /= norm2(); return *this; }
	// Normalizes the vector
	FL3 toUnitVector() const { return *this / norm2(); }
	// Any NaN?
	bool anynan() const { return isnan(x) || isnan(y) || isnan(z); }
};

class UI3
{
public:
	// Components
	size_t a, b, c;
};

#include <complex>

using FLC = std::complex<FLP>;

#include <algorithm>
#include <cctype>
#include <locale>

// trim from start (in place)
static inline void ltrim(std::string& s) {
	s.erase(s.begin(),
		std::find_if(s.begin(), s.end(),
		[](unsigned char ch) {
			return !std::isspace(ch);
		})
	);
}

// trim from end (in place)
static inline void rtrim(std::string& s) {
	s.erase(std::find_if(s.rbegin(), s.rend(),
		[](unsigned char ch) {
			return !std::isspace(ch);
		}
	).base(),s.end());
}

// trim from both ends (in place)
static inline void trim(std::string& s) {
	ltrim(s);
	rtrim(s);
}

// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
	ltrim(s);
	return s;
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
	rtrim(s);
	return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
	trim(s);
	return s;
}

// Ref: https://stackoverflow.com/a/3824338
template <class T>
void endswap(T *objp)
{
  unsigned char *memp = reinterpret_cast<unsigned char*>(objp);
  std::reverse(memp, memp + sizeof(T));
}

#endif // !__UTILITY_HPP__
