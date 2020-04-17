
#include "spherical.h"

#include <assert.h>
#include <stdlib.h>

double s2distance(double a1, double z1, double f2, double z2)
{
	double c = cos(z1) * cos(z2) + sin(z1) * sin(z2) * cos(a1 - f2);
	return acos(c);
}

void cross(vector *w, vector u, vector v)
{
	w->x = u.y * v.z - u.z * v.y;
	w->y = u.z * v.x - u.x * v.z;
	w->z = u.x * v.y - u.y * v.x;
}

void s2trinormal(vector *n, point a, point b, point c)
{
	// Convert point a to vector BA
	a.x -= b.x; a.y -= b.y; a.z -= b.z;
	// Convert point c to vector BC
	c.x -= b.x; c.y -= b.y; c.z -= b.z;
	// calculate
	cross(n, a, c);
}

void s2tricheck(s2trimesh m)
{
	// Check: V, E, F, memory allocation
	assert(m.nv > 0);
	assert(m.ne > 0);
	assert(m.nf > 0);
	assert(m.pv);
	assert(m.pe);
	assert(m.pf);

	// Validate faces
	vector normal;
	unsigned vertices[6], edge, reverse;
	for (unsigned i = 0; i < 3 * m.nf; ++i)
	{
		edge = m.pf[i] / 2;
		assert(edge < m.ne);
		reverse = !!(m.pf[i] % 2);
		vertices[(2 * i) % 6] = m.pe[edge * 2 + reverse];
		vertices[(2 * i + 1) % 6] = m.pe[edge * 2 + !reverse];

		// Vertices for a face fully loaded
		if (i % 3 == 2)
		{
			// Check: Valid edge
			assert(vertices[0] != vertices[1]);
			assert(vertices[2] != vertices[3]);
			assert(vertices[4] != vertices[5]);

			// Check: Valid cycle
			assert(vertices[1] == vertices[2]);
			assert(vertices[3] == vertices[4]);
			assert(vertices[5] == vertices[0]);
		}
	}
}

void s2trifree(s2trimesh *pm)
{
	if (pm->pv)
	{
		free(pm->pv);
		pm->pv = NULL;
		pm->nv = 0;
	}
	if (pm->pe)
	{
		free(pm->pe);
		pm->pe = NULL;
		pm->ne = 0;
	}
	if (pm->pf)
	{
		free(pm->pf);
		pm->pf = NULL;
		pm->nf = 0;
	}
}

void s2octahedron(s2trimesh *pm, unsigned level)
{
	point *pv;
	unsigned *pe, *pf;

	// Calculate number of vertices, edges and triangles
	pm->nv = (4<<(2*level)) + 2;
	pm->ne = 12<<(2*level);
	pm->nf = 8<<(2*level);

	// Preallocate memory needed
	pm->pv = (point*)malloc(pm->nv * sizeof(point));
	pm->pe = (unsigned*)malloc(2 * pm->ne * sizeof(unsigned));
	pm->pf = (unsigned*)malloc(3 * pm->nf * sizeof(unsigned));
	
	// Construct level 1
	// Vertices
	pv = pm->pv;
	s2pinit(pv++,  0,  0,  1);
	s2pinit(pv++,  1,  0,  0);
	s2pinit(pv++,  0,  1,  0);
	s2pinit(pv++, -1,  0,  0);
	s2pinit(pv++,  0, -1,  0);
	s2pinit(pv++,  0,  0, -1);
	// Edges
	pe = pm->pe;
	*(pe++) = 0; *(pe++) = 1;
	*(pe++) = 0; *(pe++) = 2;
	*(pe++) = 0; *(pe++) = 3;
	*(pe++) = 0; *(pe++) = 4;
	*(pe++) = 1; *(pe++) = 2;
	*(pe++) = 2; *(pe++) = 3;
	*(pe++) = 3; *(pe++) = 4;
	*(pe++) = 4; *(pe++) = 1;
	*(pe++) = 1; *(pe++) = 5;
	*(pe++) = 2; *(pe++) = 5;
	*(pe++) = 3; *(pe++) = 5;
	*(pe++) = 4; *(pe++) = 5;
	// Faces
	pf = pm->pf;
	*(pf++) =  0; *(pf++) =  8; *(pf++) =  3;
	*(pf++) =  2; *(pf++) = 10; *(pf++) =  5;
	*(pf++) =  4; *(pf++) = 12; *(pf++) =  7;
	*(pf++) =  6; *(pf++) = 14; *(pf++) =  1;
	*(pf++) = 16; *(pf++) = 19; *(pf++) =  9;
	*(pf++) = 18; *(pf++) = 21; *(pf++) = 11;
	*(pf++) = 20; *(pf++) = 23; *(pf++) = 13;
	*(pf++) = 22; *(pf++) = 17; *(pf++) = 15;

	// Repeatedly Refine
	for (unsigned r = 0; r < level; ++r)
	{
		assert(level == 0);
	}
}
