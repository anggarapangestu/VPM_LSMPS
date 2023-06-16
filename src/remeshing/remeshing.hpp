#ifndef INCLUDED_REMESHING
#define INCLUDED_REMESHING
#include <iomanip> // std::setw, std::setfill
#include <sstream> // std::stringstream
#include <fstream> // std::c_str()

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

#ifndef INCLUDED_LSMPS
#include "../LSMPS/LSMPSb.hpp"
#endif

#ifndef INCLUDED_NEIGHBOR
#include "../neighbor/neighbor.hpp"
#endif

class remeshing
{
	// * Internal variables
	// base_remeshing d_base_remeshing; // The method for intersearch
	std::vector<bool> sign_par;      // sign particles. @param lifetime: singleton
	Particle particleBase;           // base particle for LSMPS calculation robustness. @param lifetime: singleton
	neighbor d_neighbor;             // neighbor evaluation and base Cell List

	// Temporary internal variable for particle adaptation
	std::vector<double> adtParSize;  // Real size of particle after adaptation
	bool adap_flag;		     	 	 // The adaptation flag
	int initial_num;

	// The private method
	void redistribute_particles(Particle &p);
	LSMPSb lsmpsb;

public:
	// The public method
	void remeshing_init(Particle &p);                                 // Remeshing initialization
	bool get_remeshing(Particle &p, const Body &b, const int nIter);  // Remeshing method
};

#endif
