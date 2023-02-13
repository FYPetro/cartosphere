
#ifndef __ND_HPP__
#define __ND_HPP__

#include "cartosphere/utility.hpp"

class FL3
{
public:
	// Components
	double x, y, z;

public:
	// Default Constructor
	FL3() : x(0), y(0), z(0) {}

	// Construct from components
	FL3(double a, double b, double c) : x(a), y(b), z(c) {}

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
	FL3 operator*(double c) const { return FL3(x * c, y * c, z * c); }

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
	FL3& operator*=(double b)
	{
		x *= b;
		y *= b;
		z *= b;
		return *this;
	}

	// Component-wise division by constant
	FL3& operator/=(double b)
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
	friend FL3 operator*(double c, const FL3& b)
	{
		return b * c;
	}

	// Component-wise division by constant
	friend FL3 operator/(const FL3& a, double b)
	{
		return FL3(a.x / b, a.y / b, a.z / b);
	}

	// Dot product
	friend double dot(const FL3& a, const FL3& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	// Cross product
	friend FL3 cross(const FL3& a, const FL3& b)
	{
		double x = a.y * b.z - a.z * b.y;
		double y = a.z * b.x - a.x * b.z;
		double z = a.x * b.y - a.y * b.x;
		return FL3(x, y, z);
	}

	// Triple product
	friend double triple(const FL3& a, const FL3& b, const FL3& c)
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
	double norm2sq() const { return pow(x, 2) + pow(y, 2) + pow(z, 2); }

	// Returns the 2-norm
	double norm2() const { return sqrt(norm2sq()); }

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

#endif //!__ND_HPP__
