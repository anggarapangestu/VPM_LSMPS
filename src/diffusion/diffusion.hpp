#ifndef INCLUDED_DIFFUSION
#define INCLUDED_DIFFUSION

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

/**
 *  @brief A particle diffusion subroutine. Perform diffusion to the particle
 *  only for the active region. Also update the vorticity value.
 *
 *  @headerfile diffusion.hpp
 */
class diffusion
{
	// A temporary method
	std::vector<double> dvordt;	// The collection of the particle vorticity 
	std::vector<double> dvordt_prev;	// The previous value on simulation
	Particle activePar;
	std::vector<int> index;

	// Diffusion calculation
	void diffusion_2d(Particle &_particle);
	void diffusion_3d(Particle &_particle);
	void identify_active_particle(Particle &_particle);

	void diffusion_2d_corr(Particle &_particlePred, Particle &_particle);
	void diffusion_3d_corr(Particle &_particlePred, Particle &_particle);		// Not yet created

public:
	// Diffusion manager calculation
	void main_diffusion(Particle &_particle);
	// void diffusion_prediction(Particle &_particle, std::vector<double>& dwdt);
	void diffusion_correction(Particle &_particlePred, Particle &_particle);
	void diffusion_AB2(Particle &_particle);
};

#endif
