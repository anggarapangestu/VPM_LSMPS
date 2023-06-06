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

	// // TODO: PSE scheme
	// // ! change input [pair_i, pair_j]
	// std::vector<double> pse(int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &sp, const std::vector<double> &gpz, const std::vector<int> &pair_i, const std::vector<int> &pair_j);

public:
	// Diffusion calculation
	void main_diffusion(Particle &p);
};

#endif
