
#ifndef __SPHERICAL_H__
#define __SPHERICAL_H__

#include <math.h>

static inline double deg2rad(double degrees)
{
	return (M_PI / 180.f) * degrees;
}

static inline double rad2deg(double radians)
{
	return (M_1_PI * 180.f) * radians;
}

static inline double zen2lat(double zenith)
{
	return M_PI_2 - zenith;
}

static inline double lat2zen(double latitude)
{
	return M_PI_2 - latitude;
}

static inline double azi2lon(double azimuth)
{
	double normalized = azimuth / (2 * M_PI) - 0.5;
	return azimuth - (2 * M_PI) * ceil(normalized);
}

double s2distance(double, double, double, double);

struct double3 {
	double x;
	double y;
	double z;
};

static inline double dot(struct double3 u, struct double3 v)
{
	return u.x * v.x + u.y * v.y + u.z * v.z;
}

static inline void s2pinit(struct double3 *p, double x, double y, double z)
{
	p->x = x;
	p->y = y;
	p->z = z;
}

static inline void s2xyz(struct double3 *p, double a, double z)
{
	double r = sin(z);
	p->z = cos(z);
	p->x = r * cos(a);
	p->y = r * sin(a);
}

typedef struct double3 point;

typedef struct double3 vector;

void cross(vector *w, vector u, vector v);

typedef struct {
	unsigned nv;
	unsigned ne;
	unsigned nf;
	point *pv;
	unsigned *pe;
	unsigned *pf;
} s2trimesh;

void s2tricheck(s2trimesh);

void s2trifree(s2trimesh *);

void s2octahedron(s2trimesh *, unsigned);

#endif //!__SPHERICAL_H__
