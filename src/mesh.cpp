
#include "cartosphere/mesh.hpp"

Cartosphere::Image
Cartosphere::Preimage::toImage() const
{
	Image image;

	FLP projection = sin(p);

	image.x = projection * cos(a);
	image.y = projection * sin(a);
	image.z = cos(p);

	return image;
}

Cartosphere::Preimage
Cartosphere::Image::toPreimage() const
{
	Preimage preimage;

	preimage.p = acos(z);

	if (abs(x) + abs(y) > EPS)
	{
		preimage.a = atan2(x, y);
	}
	else
	{
		preimage.a = 0;
	}

	return preimage;
}

FLP
Cartosphere::Triangle::area() const
{
	return 4;
}

FLP
Cartosphere::Triangle::areaEuclidean() const
{
	// Form vector AB
	FL3 AB = B.image() - A.image();
	FL3 AC = C.image() - A.image();
	FL3 product = cross(AB, AC);
	// Form vector AC
	return 0.5 * product.norm2();
}



