#ifndef PARTICLE_ADAPTATION
#define PARTICLE_ADAPTATION

#include "../grid_block/gridNodeAdapt.hpp"

class adaptation
{
private:
    // Internal variables
    std::vector<std::vector<double>*> predProp;      // The container of properties to be considered
    
    // * Internal variables
    void error_selector(const Particle &_parEval);      // Select the error predictor type (selection on Pars namespace)
    void error_adjustment(int winPos);                            // Adjust the error selector
    void feature_predictor(const Particle &_parEval);   // Base on the 0th order of vorticity value
    void gradient_predictor(const Particle &_parEval);  // Base on the 2nd order vorticity gradient

    void laplace_predictor(const Particle &_parEval);  // Base on the 2nd order vorticity gradient
public:
    bool get_adaptation(Particle &_parEval, Particle *&_parBase, GridNode &_baseGrid);

    void get_adaptation_LSMPS(Particle &_parEval, GridNode &_baseGrid);
};

#endif
