#include "diffusion.hpp"

void diffusion::main_diffusion(Particle &p)
{
    // Internal variables
    Particle _particle;           // The data of particle inside 'active particle + buffer zone'
    std::vector<double> dvordt;   // The vorticity time differential
    std::vector<int> _index;      // The index list of particle inside 'active particle + buffer zone'
    std::vector<bool> _eval;      // The flag list of evaluated particle
    _eval = std::vector<bool>(p.num, false);    // At initial still no evaluated particle list
    
    // Diffusion computational time intialization
    clock_t t, _t;
    t = clock();        // Total computation time
    _t = clock();       // Local computation time
    
    // Diffusion prompt
    printf("\nCalculating diffusion ...\n");
    
    // Assign the index of active particles
    for (size_t ID = 0; ID < p.num; ID++)
    {
        // Check if the particle ID is active
        if (p.isActive[ID] == true)
        {
            // [1] The current ID particle
            if(_eval[ID] == false)
            {
                // Assign ID into index list if not in the list (if flag == false)
                _index.push_back(ID);    // Assign ID into the index list
                _eval[ID] = true;        // Update the evaluation flag of particle ID
            }

            // [2] The neighbor particles of current ID particle
            for (auto _ngh_ID : p.neighbor[ID])
            {
                if(_eval[_ngh_ID] == false)
                {
                    // Assign into index list if not in the list (if flag == false)
                    _index.push_back(_ngh_ID);  // Assign _ngh_ID particle into the index list
                    _eval[_ngh_ID] = true;      // Update the evaluation flag of particle _ngh_ID
                }
            }
        }
    }

    // // [OLD ALGORITHM] Store only unique ID 
    // std::sort(_index.begin(), _index.end());
    // std::vector<int>::iterator it;
    // it = std::unique(_index.begin(), _index.end());
    // _index.resize(std::distance(_index.begin(), it));

    // Display local computational time
    _t = clock() - _t;
    printf("<-> Evaluating partilce for diffusion: [%f s]\n", (double)_t/CLOCKS_PER_SEC);
    _t = clock();

    // Resize the _particle
    _particle.num = _index.size();
    _particle.x.resize(_particle.num);
    _particle.y.resize(_particle.num);
    _particle.s.resize(_particle.num);
    _particle.vorticity.resize(_particle.num);
    _particle.neighbor.resize(_particle.num);
    
    // Store the particle data of each particle inside _index list
    for (size_t i = 0; i < _particle.num; i++)
    {
        int _ID = _index[i];
        _particle.x[i] = (p.x[_ID]);
        _particle.y[i] = (p.y[_ID]);
        _particle.s[i] = (p.s[_ID]);
        _particle.vorticity[i] = (p.vorticity[_ID]);
        _particle.neighbor[i] = (p.neighbor[_ID]);
    }

    // Calculate the second order differential of vorticity
    lsmpsa.set_LSMPS(_particle.x, _particle.y, _particle.s, _particle.vorticity, p.x, p.y, p.s, p.vorticity, _particle.neighbor);
    std::vector<double> _d2fd2x = lsmpsa.get_d2d2x();
    std::vector<double> _d2fd2y = lsmpsa.get_d2d2y();
    
    // Calculate the vorticity laplacian to get the vorticity time differential
    dvordt.resize(p.num, 0.0e0);
    for (int i = 0; i < _particle.num; i++)
    {
        dvordt[_index[i]] = Pars::NU * ((_d2fd2x[i] + _d2fd2y[i])); 
    }

    // Diffusion time integration by 1st order; [HYPOTHESIS: consider the diffusion only]
    for (size_t i = 0; i < p.num; i++)
    {
        p.vorticity[i] += Pars::dt * (dvordt[i]);                // w = del x u
        p.gz[i] += Pars::dt * (dvordt[i]) * std::pow(p.s[i], 2); // γ = w * Area
    }

    // Display local computational time
    _t = clock() - _t;
    printf("<-> Calculating diffusion:             [%f s]\n", (double)_t/CLOCKS_PER_SEC);

    // Display computational time
    t = clock() - t;
    printf("<-> Diffusion total comp. time:        [%f s]\n", (double)t/CLOCKS_PER_SEC);
}
