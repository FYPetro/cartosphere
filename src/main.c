
#include <stdio.h>
#include <stdlib.h>

#include "spherical.h"

int main(int argc, char *argv[])
{
	s2trimesh mesh;

	printf("INITIALIZING...\n");
	s2octahedron(&mesh,0);
	s2tricheck(mesh);

	printf("FREEING...\n");
	s2trifree(&mesh);

	printf("DONE...\n");
	return 0;
}
