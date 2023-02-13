
#include "cartosphere/functions.hpp"

#ifdef APPLE_LIKE
#include <boost/math/special_functions/spheric_harmonic.hpp>
using boost::math::spherical_harmonic_r;
using boost::math::spherical_harmonic_i;
#endif

double
cs_y(int l, int m, double theta, double phi)
{
	double r;
#ifdef APPLE_LIKE
	if (m >= 0)
	{
		r = spherical_harmonic_r(l, m, theta, phi);
	}
	else
	{
		r = spherical_harmonic_i(l, -m, theta, phi);
	}
	if (m != 0)
	{
		r *= M_SQRT2;
	}
#else
	r = std::sph_legendre(l, abs(m), theta);
	if (m != 0)
	{
		// Remove the Condon-Shortley phase and provide real normlization
		r *= (m % 2) ? (-M_SQRT2) : M_SQRT2;
		// Attach the azimuthal factor
		if (m > 0)
		{
			r *= cos(m * phi);
		}
		else
		{
			r *= sin((-m) * phi);
		}
	}
#endif
	return r;
}

double
cs_azview(double p, double a, double pv, double av)
{
	double v, t, b;

	if (pv == 0)
	{
		v = a;
	}
	else
	{
		t = sin(av - a) * sin(p);
		b = cos(p) * sin(pv) - sin(p) * cos(pv) * cos(av - a);
		v = atan2(t, b);
	}

	return v;
}
