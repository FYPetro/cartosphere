
#ifndef __GLOBE_HPP__
#define __GLOBE_HPP__

#include "cartosphere/mesh.hpp"

#include <string>

namespace Cartosphere
{

	class Globe
	{
	public:
		// Construct a globe from file
		Globe(const std::string& path);

	public:
		// Timestepping
		void run();

		// Output the result
		void format(const std::string& path);

	private:
		// Internal computational mesh
		Cartosphere::TriangularMesh _m;
	};

}

#endif // !__GLOBE_HPP__
