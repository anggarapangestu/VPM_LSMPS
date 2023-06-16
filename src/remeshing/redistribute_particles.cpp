#include "remeshing.hpp"
#include <omp.h>

// Redistribution particle using LSMPS B interpolation on vorticity
void remeshing::redistribute_particles(Particle &parEval)
{   
    // Here the procedure to add the current particle as the neighbor need to be discussed after the particle adaptation is work
    clock_t t;
    t = clock();
    
    // TODO: Assign the neighbor
    std::vector<std::vector<int>> _par2grdNeighbor;  // The neighbor of base particle using the ID of new particle (after advected)
    std::vector<std::vector<int>> _grd2grdNeighbor;  // The neighbor of base particle relative to itself
    
    // Evaluate the _grd2grd neighbor list
    _grd2grdNeighbor = this->particleBase.neighbor;
    // Add the particle ID itself into the neighborID
    for (int i = 0; i < this->initial_num; i++){
        _grd2grdNeighbor[i].push_back(i);
    }

    // Evaluate the _par2grd neighbor list
    bool parAdapt = false;
    if (parAdapt == false){
        _par2grdNeighbor = this->particleBase.neighbor;
        // Add the particle ID itself into the neighborID
        for (int i = 0; i < this->initial_num; i++){
            _par2grdNeighbor[i].push_back(i);
        }
    }

	t = clock() - t;
    printf("<-> Evaluating partilce neighbor list: [%f s]\n", (double)t/CLOCKS_PER_SEC);
    t = clock();
    
    // Update the _par2grdNeighbor based on the first neighbor particle
    // this->d_neighbor.inter_search(parEval, this->particleBase, _par2grdNeighbor);
    // _par2grdNeighbor = this->d_neighbor.par2grid_neigbor_search(parEval, this->particleBase);

    // TODO: redistribute particles
    // Perform LSMPS interpolation of vorticity
    lsmpsb.set_LSMPS(this->particleBase.x, this->particleBase.y, 
                     this->particleBase.s, this->particleBase.vorticity, 
                     parEval.x, parEval.y, parEval.s, parEval.vorticity,
                     _par2grdNeighbor/*_grd2grdNeighbor*/, _grd2grdNeighbor);
    
    t = clock() - t;
    printf("<-> Calculating LSMPS:                 [%f s]\n", (double)t/CLOCKS_PER_SEC);
    t = clock();

    // Update the particle base data vorticity and vortex strength
    this->particleBase.vorticity = lsmpsb.get_d00();
    if (this->adap_flag == true){
        this->particleBase.s = this->adtParSize;
    }
    for (size_t _i = 0; _i < this->particleBase.num; _i++){
        this->particleBase.gz[_i] = this->particleBase.vorticity[_i] * std::pow(this->particleBase.s[_i],2);
    }

    // TODO: Update active particle
    double vor_max = 0.0e0;
    for (int i = 0; i < this->particleBase.num; i++)
    {
        vor_max = vor_max > std::abs(this->particleBase.vorticity[i]) ? 
                  vor_max : std::abs(this->particleBase.vorticity[i]);
    }

    for (size_t i = 0; i < this->particleBase.num; i++)
    {
        double vor = std::abs(this->particleBase.vorticity[i]);
        if (vor >= 1.0e-4 * vor_max || this->particleBase.isActive[i] == true
            // vor >= Pars::re_trsh * (Pars::vis * p.s[i]) ||
            )
        {
            // Only consider active if it is not 
            // if (this->particleBase.chi[i] != 1.0){
                this->particleBase.isActive[i] = true;
            // }
        }
        // else if (vor >= 1.0e-7 * vor_max){
        //     // this->particleBase.isActive[i] = true;
        // }
        else
        {
            this->particleBase.isActive[i] = false;
            // prevent truncated error
            this->particleBase.vorticity[i] = 0.0e0;
            this->particleBase.gz[i] = 0.0e0;  
        }
    }

    // Update the redistribute data into the particle variable
    parEval = this->particleBase;
    // parEval.num = this->particleBase.num;
    // parEval.neighbor = this->particleBase.neighbor;
    // parEval.basis_label = this->particleBase.basis_label;
    // parEval.cell_label = this->particleBase.cell_label;
    // parEval.x = this->particleBase.x;
    // parEval.y = this->particleBase.y;
    // parEval.s = this->particleBase.s;
    // parEval.gz = this->particleBase.gz;
    // parEval.vorticity = this->particleBase.vorticity;
    // parEval.isActive = this->particleBase.isActive;


    t = clock() - t;
    printf("<-> Re-assigning the particle data:    [%f s]\n", (double)t/CLOCKS_PER_SEC);
}

/* [OLD CODE for brute force enlarging neighbor particle]
for (int i = 0; i < this->particleBase.num; i++){
    std::vector<int> _tempList;
    for (auto ID:parEval.neighbor[i]){
        // List of all neighbor of the current particle
        for (auto _ID:parEval.neighbor[ID]){
            // List of all neighbor ID of each neighbor particle
            _tempList.emplace_back(_ID);
        }
    }
    
    // Store only unique ID
    std::sort(_tempList.begin(), _tempList.end());
    std::vector<int>::iterator it;
    it = std::unique(_tempList.begin(), _tempList.end());
    _tempList.resize(std::distance(_tempList.begin(), it));
    
    // Assign the new neighbor ID
    _par2grdNeighbor.emplace_back(_tempList);
}
*/