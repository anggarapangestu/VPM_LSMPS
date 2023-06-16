#ifndef INCLUDED_DIFFUSION
#define INCLUDED_DIFFUSION

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

// #ifndef INCLUDED_BIOT_SAVART
// #include "../velocity_poisson/velocity_biot_savart.hpp"
// #endif

#ifndef INCLUDED_LSMPSa
#include "../LSMPS/LSMPSa.hpp"
#endif

// class velocity_biot_savart;
// class base_dc;
// class neighbor;
class LSMPSa;
class diffusion
{
	// Internal variable
	// velocity_biot_savart d_base_poisson;
	LSMPSa lsmpsa;		 // To calculate laplacian

public:
	// Diffusion calculation
	void main_diffusion(Particle &p);
};

#endif
