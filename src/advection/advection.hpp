#ifndef INCLUDED_ADVECTION
#define INCLUDED_ADVECTION

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

#include "../grid_block/gridNodeNgh.hpp"		// Grid node method

/**
 *  @brief A particle advection subroutine. Perform advection to the particle
 *  using the given type. There are euler and runge kutta 2 type.
 *
 *  @headerfile advection.hpp
 */
class advection
{
	// The euler advection
	void advection_euler(Particle &_particle);
	
	// Trying using runge kutta 2nd order
	void advection_rk2(Particle &_particle);
public:
	// The public advection
	void main_advection(Particle &_particle);
	
	void main_advection_corr(Particle &_particlePred, Particle &_particle);

	// The advection calculation using the Adam Bashfort Time Integration
	void main_advection_AB(Particle &_particlePrev, Particle &_particle);

	// Interpolate the velocity into new particle location
	void intplt_velocity(Particle &_particleInit, Particle &_particleAdv, const GridNode &_baseGrid);
};

#endif
