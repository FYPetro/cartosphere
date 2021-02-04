
#ifndef __UTILITY_HPP__
#define __UTILITY_HPP__

using FLP = double;

template<typename T>
constexpr T deg2rad(T angle) { return (angle * M_PI / (FLP)180); }

template<typename T>
constexpr T rad2deg(T angle) { return (angle * M_1_PI * (FLP)180); }

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
	FL3(FL3 const& that) = default;
	// Assignment operator
	FL3& operator=(FL3 const& that)
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
	// Component-wise division by constant
	FL3& operator/=(FLP b)
	{
		x /= b;
		y /= b;
		z /= b;
		return *this;
	}
	// Component-wise addition
	friend FL3 operator+(FL3 const& a, FL3 const& b)
	{
		return FL3(a.x + b.x, a.y + b.y, a.z + b.z);
	}
	// Component-wise subtraction
	friend FL3 operator-(FL3 const& a, FL3 const& b)
	{
		return FL3(a.x - b.x, a.y - b.y, a.z - b.z);
	}
	// Scalar multiplication
	friend FL3 operator*(FLP c, FL3 const& b)
	{
		return b * c;
	}
	// Component-wise division by constant
	friend FL3 operator/(FL3 const& a, FLP b)
	{
		return FL3(a.x / b, a.y / b, a.z / b);
	}
	// Cross product
	friend FL3 cross(FL3 const& a, FL3 const& b)
	{
		FLP x = a.y * b.z - a.z * b.y;
		FLP y = a.z * b.x - a.x * b.z;
		FLP z = a.x * b.y - a.y * b.x;
		return FL3(x, y, z);
	}
	// Normalize
	friend FL3 normalize(FL3 const& a)
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
};

class UI3
{
public:
	// Components
	size_t a, b, c;
};

#endif // !__UTILITY_HPP__
