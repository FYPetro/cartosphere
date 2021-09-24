
#include "cartosphere/functions.hpp"

#include <cmath>

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
