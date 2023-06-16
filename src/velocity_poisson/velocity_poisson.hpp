#ifndef INCLUDED_VELOCITY_POISSON
#define INCLUDED_VELOCITY_POISSON

#include <string>

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

#ifndef INCLUDED_BIOT_SAVART
#include "velocity_biot_savart.hpp"
#endif

#ifndef INCLUDED_LSMPSa
#include "../LSMPS/LSMPSa.hpp"
#endif

#ifndef INCLUDED_LSMPSb
#include "../LSMPS/LSMPSb.hpp"
#endif

#ifndef INCLUDED_NEIGHBOR
#include "../neighbor/neighbor.hpp"
#endif

#ifndef FMM_2D_PACKAGE
#include "../FMM2D/fmm2D.hpp"
#endif

class LSMPSa;
class LSMPSb;
class neighbor;
class velocity_biot_savart;
class base_remeshing;

class velocity_poisson
{
private:
	// -- creating instance
	velocity_biot_savart biot_savart;
	neighbor d_neighbor;			// linked_list, direct_find
	// base_remeshing _base_remeshing; // inter-particle searching
	LSMPSb _lsmpsb;
	treeCell treeData;				// The basis of tree data of poisson solver

public:
	void get_velocity(Particle &p, const int t); //udah diedit, awalnya void poisson (Particle &p)
};

#endif
