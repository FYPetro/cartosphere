
#ifndef __UTILITY_HPP__
#define __UTILITY_HPP__

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>

#define FLP double

static const FLP EPS = std::numeric_limits<FLP>::epsilon();

// double deg2rad(double degree)
// {
// 	return degree * M_PI / 180;
// }

// double rad2deg(double radian)
// {
// 	return radian / M_PI * 180;
// }

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
	FL3(FL3 const &that) = default;
	// Assignment operator
	FL3 &operator=(FL3 const &that)
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
	// Vector subtraction
	friend FL3 operator-(FL3 const &a, FL3 const &b)
	{
		return FL3(a.x - b.x, a.y - b.y, a.z - b.z);
	}
	// Cross product
	friend FL3 cross(FL3 const &a, FL3 const &b)
	{
		FLP x = a.y * b.z - a.z * b.y;
		FLP y = a.z * b.x - a.x * b.z;
		FLP z = a.x * b.y - a.y * b.x;
		return FL3(x, y, z);
	}

public:
	// Returns the square of the 2-norm
	FLP norm2sq() const { return pow(x,2) + pow(y,2) + pow(z,2); }
	// Returns the 2-norm
	FLP norm2() const { return sqrt(norm2sq()); }
};

class UI3
{
public:
	// Components
	unsigned a, b, c;
};

#endif // !__UTILITY_HPP__
