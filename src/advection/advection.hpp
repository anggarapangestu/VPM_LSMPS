#ifndef INCLUDED_ADVECTION
#define INCLUDED_ADVECTION

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

class advection
{
public:
	// The basic advection
	void advection_euler(Particle &p);
	
	// Trying using runge kutta 2nd order
	void advection_rk2(Particle &p, std::vector<double> &dfdt);
	
	// void advection(Particle &p, const Particle &pBase, std::vector<double> &dfdt);
};

#endif
