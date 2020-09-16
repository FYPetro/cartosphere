
#ifndef _PREDICATES_HPP
#define _PREDICATES_HPP

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif //!_USE_MATH_DEFINES

#include <cmath>
#include <limits>

#define FLP double

// List of spherical geometry predicates
#ifndef _DISABLE_PREDICATES

template<typename T>
constexpr T deg2rad(T angle) { return (angle * M_PI / (FLP)180); }

template<typename T>
constexpr T rad2deg(T angle) { return (angle * M_1_PI * (FLP)180); }

static const FLP EPS = std::numeric_limits<FLP>::epsilon();

namespace Cartosphere
{
	namespace Predicates
	{
		template<typename T>
		constexpr auto is_pi(T angle) { return angle == M_PI; }
	}
}
#endif //!_DISABLE_PREDICATES

#endif //!_PREDICATES_HPP
