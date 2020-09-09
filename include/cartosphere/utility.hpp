
#ifndef __UTILITY_HPP__
#define __UTILITY_HPP__

#ifndef _DISABLE_PREDICATES
#define _DISABLE_PREDICATES
#include "cartosphere/predicates.hpp"
#undef _DISABLE_PREDICATES
#endif //!_DISABLE_PREDICATES

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
	// Component-wise division by constant
	FL3 &operator/=(FLP b)
	{
		x /= b;
		y /= b;
		z /= b;
		return *this;
	}
	// Component-wise addition
	friend FL3 operator+(FL3 const &a, FL3 const &b)
	{
		return FL3(a.x + b.x, a.y + b.y, a.z + b.z);
	}
	// Component-wise subtraction
	friend FL3 operator-(FL3 const &a, FL3 const &b)
	{
		return FL3(a.x - b.x, a.y - b.y, a.z - b.z);
	}
	// Component-wise division by constant
	friend FL3 operator/(FL3 const &a, FLP b)
	{
		return FL3(a.x / b, a.y / b, a.z / b);
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
	// Normalizes the vector
	FL3 &normalize() { (*this) /= norm2(); return *this; }
};

class UI3
{
public:
	// Components
	unsigned a, b, c;
};

#endif // !__UTILITY_HPP__
