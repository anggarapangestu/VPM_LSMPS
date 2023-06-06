#ifndef INCLUDED_PENALIZATION
#define INCLUDED_PENALIZATION
#include <fstream> // data stream

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

#ifndef INCLUDED_GEOMETRY
#include "../geometry/geometry.hpp"
#endif

#ifndef INCLUDED_NEIGHBOR
#include "../neighbor/neighbor.hpp"
#endif

#ifndef INCLUDED_SAVE_DATA_BASE
#include "../save_data/base_save_data.hpp"
#endif

#ifndef INCLUDED_LSMPSa
#include "../LSMPS/LSMPSa.hpp"
#endif

// #ifndef INCLUDED_FORTRAN_UTILS
// #include "../Fortran/Utils/FortranUtils.hpp"
// #endif

class penalization
{
	// ** Internal variables
	double lambda;			 // porosity parameter

	// ** Creating instance
	neighbor d_neighbor;
	geometry d_geom;
	base_save_data d_base_save_data;

	// Internal Method
	void kai_def(const Body &b,  Particle &_p);  // Calculate kai parameter
	void lambda_def();                           // Define the value of lambda
	void no_slip(Particle &_p, const Body &b);

	// Additional
	void no_slip_iterative(Particle &_p, Particle &p, const Body &b, int it);

public:
	// Main penalization manager
	void get_penalization(Particle &p, const Body &b, int it);
};

#endif
