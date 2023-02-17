
#ifndef __FUNCTIONS_HPP__
#define __FUNCTIONS_HPP__

// Use Cartosphere's own utility tools
#include "cartosphere/utility.hpp"

// Use core math functionality
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !_USE_MATH_DEFINES
#include <cmath>

// Character types to support string manipulations
#include <cctype>
#include <locale>

// Defines associated Legendre polynomial of degree l and order m without the
// (odd) Condon-Shortley phase factor
// [C++17] Using fold expressions to pass along arguments
// Ref: https://en.cppreference.com/w/cpp/language/fold
#ifdef APPLE_LIKE
#include <boost/math/special_functions/legendre.hpp>
template <typename T>
constexpr T cs_legendre(int l, T x)
{
	return boost::math::legendre_p(l, x);
}

template <typename T>
constexpr T cs_legendre(int l, int m, T x)
{
	return ((m % 2) ? -1 : 1) * boost::math::legendre_p(l, m, x);
}
#else
constexpr auto cs_legendre = [](auto &&...args) {
	return std::legendre(std::forward<decltype(args)>(args)...);
};
#endif

// Orthonormal real-valued spherical harmonic basis function
double cs_y(int l, int m, double theta, double phi);

// Calculate the azimuth of one point from the point of view of another
double cs_azview(double p, double a, double pv, double av);

// Convert degrees to radians
template <typename T>
constexpr double cs_deg2rad(T angle) { return (angle * M_PI / 180); }

// Convert radians to degrees
template <typename T>
constexpr double cs_rad2deg(T angle) { return (angle * M_1_PI * 180); }

// Trim from start (in place)
static inline void ltrim(string& s) {
	s.erase(s.begin(),
		std::find_if(s.begin(), s.end(),
			[](unsigned char ch) {
				return !std::isspace(ch);
			})
	);
}

// Trim from end (in place)
static inline void rtrim(string& s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(),
		[](unsigned char ch) {
			return !std::isspace(ch);
		}
	).base(), s.end());
}

// Trim from both ends (in place)
static inline void trim(string& s)
{
	ltrim(s);
	rtrim(s);
}

// Trim from start (copying)
static inline string ltrim_copy(string s)
{
	ltrim(s);
	return s;
}

// Trim from end (copying)
static inline string rtrim_copy(string s)
{
	rtrim(s);
	return s;
}

// Trim from both ends (copying)
static inline string trim_copy(string s)
{
	trim(s);
	return s;
}

// Switch endian
// By default, little endian is used on Windows
// Ref: https://stackoverflow.com/a/3824338
template <typename T>
void endswap(T* objp)
{
	unsigned char* memp = reinterpret_cast<unsigned char*>(objp);
	std::reverse(memp, memp + sizeof(T));
}

#endif // !__FUNCTIONS_HPP__
