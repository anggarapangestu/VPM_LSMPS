#ifndef INCLUDED_STR
#define INCLUDED_STR

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

/**
 *  @brief A particle stretching subroutine. Perform the stretching to the particle
 *  only for the active region. Also update the vorticity value.
 *
 *  @headerfile stretching.hpp
 */
class stretching
{
    // A temporary method
	std::vector<double> dvorXdt;
    std::vector<double> dvorYdt;
    std::vector<double> dvorZdt;

    std::vector<double> dvorXdtPrev;
    std::vector<double> dvorYdtPrev;
    std::vector<double> dvorZdtPrev;

	Particle activePar;
	std::vector<int> index;

    void identify_active_particle(Particle &_particle);

public:
    // Stretching calculation
	void main_stretching(Particle &_particle);

    void main_stretching_corr(Particle &_particlePred, Particle &_particle);

    // void main_stretching_AB(Particle &_particle, std::vector<std::vector<double>*> &VortDiff);
    void main_stretching_AB(Particle &_particle);

    // Combination of diffusion and stretching only for 3D simulation
	void calc_diff_stretch(Particle &_particle);

    // Combination of diffusion and stretching only for 3D simulation (Using to initiate the Adam-Bashfort)
    void calc_diff_stretch(Particle &_particle, std::vector<std::vector<double>*> &VortDiff);       // The combination of adam-bashfort and the predictor corrector (Actually Forward Euler)

    // Combination of diffusion and stretching only for 3D simulation for Adam-Bashfort
    void calc_diff_stretch_AB(Particle &_particle, std::vector<std::vector<double>*> &VortDiff);

    void diff_stretch_correction(Particle &_particlePred, Particle &_particle);

};

#endif
