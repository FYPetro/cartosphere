
#include "cartosphere/functions.hpp"

#include <cmath>

#ifndef _WIN32

FLP
cartosphere_Y_real(int l, int m, FLP z, FLP a)
{
	return 0;
}

FLP
cartosphere_azview(FLP p, FLP a, FLP pv, FLP av)
{
	return 0;
}

#elif _WIN32

FLP
cartosphere_Y_real(int l, int m, FLP z, FLP a)
{
	FLP r = std::sph_legendre(l, std::abs(m), z);

	if (m != 0)
	{
		r *= std::sqrt(2) * std::pow(-1, m);

		if (m > 0)
		{
			r *= cos(m * a);
		}
		else
		{
			r *= sin(m * a);
		}
	}

	return r;
}

FLP
cartosphere_azview(FLP p, FLP a, FLP pv, FLP av)
{
	FLP v, t, b;

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

#else

#endif