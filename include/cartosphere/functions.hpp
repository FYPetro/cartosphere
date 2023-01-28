
#ifndef __FUNCTIONS_HPP__
#define __FUNCTIONS_HPP__

#include "cartosphere/utility.hpp"

FLP cartosphere_Y_real(int l, int m, FLP z, FLP a);

FLP cartosphere_azview(FLP p, FLP a, FLP pv, FLP av);

// Mathematical utility
template<typename T>
constexpr FLP deg2rad(T angle) { return (angle * (FLP)M_PI / (FLP)180); }

template<typename T>
constexpr FLP rad2deg(T angle) { return (angle * (FLP)M_1_PI * (FLP)180); }

#endif // !__FUNCTIONS_HPP__
