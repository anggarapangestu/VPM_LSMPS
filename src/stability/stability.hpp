#ifndef INCLUDED_STABILITY
#define INCLUDED_STABILITY

#include "../../Utils.hpp"

/**
 *  @brief A property set of single body made of node or panel.
 *
 *  @headerfile Utils.hpp
 */
class stability
{
private:
    // The label of container consisted in 3 value given below
    // [(1)Current Evaluated -> Average value, (2)Accumulative in Domain, (3)Global Maximum]

    /** The criteria for stability to control
     *  -------------------------------------
     *    Courant Number     :     C = u * dt / h          < 1.0
     *    Diffusion Number   :   Phi = nu * dt / h^2       < 0.5
     *    Vortex Number      :  Re_h = w * h^2 / nu      [in O(1)]
     *                                  or
     *                          Re_h = w * dt            [in O(1)] (combined with diffusion)
     *    Lagrangian Stretch :   C_L = |del(u)|_max * dt   < 0.5
    */
    
    void courant_eval(const Particle _par, double* _container);          	// Evaluation on the courant number (Must below 1.0)
    void diffusion_eval(const Particle _par, double* _container);      	 	// Evaluation on the diffusion number (Must below 0.25 or 0.5) Actually diffusion is already glabaly calculated for single resolution
    void vortex_eval(const Particle _par, double* _container, int _type); 	// Evaluation on the vortex number (Must in order of O(1))
    void lag_str_eval(const Particle _par, double* _container);          	// Evaluation on the lagrangian stretching criteria (Must below value of 0.5)

public:
    void stabilityEval(const Particle &_particle, 
					   std::vector<double*> &_maxValue, 
					   const int &step);
    // The evaluation of stability consisted of maximum value, average calculation, and time simulation value
                       
};

#endif