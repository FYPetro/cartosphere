
#ifndef __GLOBE_HPP__
#define __GLOBE_HPP__

#include "cartosphere/mesh.hpp"

namespace Cartosphere
{
	class Globe
	{
	public:
		// Construct a globe from file
		Globe(const string& path);

	public:
		// Timestepping
		void run();

		// Output the result
		void format(const string& path);

	private:
		// Internal computational mesh
		Cartosphere::TriangularMesh _m;
	};
}

#endif // !__GLOBE_HPP__
