#include "neighbor.hpp"

// The particle neighbor list
void neighbor::cell_list_init(Particle& parEval){
    // Generate the Cell List
    printf("Generating cell list ...\n");
    clock_t t = clock();

    // Initialize the Cell List
    d_cell_list.initCellList(parEval);
    
    // Evaluate the particle cell ID
    d_cell_list.createCellList(parEval);

    // // Show particle neighbor
    // d_cell_list.showBasisNgh();
    
    // particle generation initialization summary time display
    t = clock() - t;
	printf("<-> Cell list initialization\n");
    printf("    comp. time:                        [%f s]\n", (double)t/CLOCKS_PER_SEC);

    // Evaluate the particle cell ID
    if (Pars::flag_save_cell){d_cell_list.saveCellData();}
    
}

// ************************************************************************************
// ====================================================================================

// The particle neighbor list
void neighbor::neigbor_search(Particle& parEval){
    // The terminology
    // * parEval is the particle to be evaluated
    //   > each of particle inside parEval will be evaluated by the parList to calculate neigbor
    //   > here is the particle where the calculated neighbor is stored
    // * parList only store the data of particle and its corresponding location
    //   > the particle inside parList store the data for neighbor position evaluation
    // * The particle of parList and parEval could be similar or different

    // --- Direct find;
    printf("Evaluating neighbor ...\n");
    clock_t t = clock();
    if (Pars::opt_neighbor == 0)
    {
        printf("<+> Type 0: Direct neighbor search\n");
        parEval.neighbor = this->direct_find(parEval.num, parEval.s, parEval.x, parEval.y, Pars::r_sup);
    }
    // --- Linked list
    else if (Pars::opt_neighbor == 1)
    {
        printf("<+> Type 1: Link list neighbor search\n");
        parEval.neighbor = this->link_list(parEval.num, parEval.s, parEval.x, parEval.y, int(Pars::r_sup));
    }
    // --- Included neighbor search package from Cell List
    else if (Pars::opt_neighbor == 2)
    {
        printf("<+> Type 2: Cell list neighbor search\n");
        this->d_cell_list.findNeighbor(parEval);
        // this->d_cell_list.checkNGH(parEval);
    }
    // --- Spatial hash
    else if (Pars::opt_neighbor == 3)
    {
        printf("<+> Type 3: Spatial hash neighbor search\n");
        this->spatial_hash(parEval);
    }
    
    // particle generation initialization summary time display
    t = clock() - t;
    printf("<-> Neighbor search comp. time:        [%f s]\n", (double)t/CLOCKS_PER_SEC);

}

// ************************************************************************************
// ====================================================================================

// Particle adaptation 
bool neighbor::particle_adaptation(const Particle& parEval, Particle& baseParticle, std::vector<double>& PARsize){
    // Note: At initial the distribution of parEval and baseParticle is the same but not the properties

    // computational time accumulation
    clock_t t;
    t = clock();
    /* Procedure of particle adaptation:
       [1] Evaluate each particle, determine if need adaptation of not
       [2] Perform the particle adaptation if necesarry (based on the procedure 1)
    */

    // PROCEDURE 1 ! : Evaluating particle adaptation
    // *************
    // The marker whether the adaptation need to performed or not
    bool _adaptation = false;

    // TODO: Find maximum vorticity
    double vor_max = 0.0e0;
    for (int i = 0; i < parEval.num; i++)
    {
        vor_max = vor_max > std::abs(parEval.vorticity[i]) ? 
                  vor_max : std::abs(parEval.vorticity[i]);
    }

    // Initialize the Cell Level Up setting (Setting for particle adaptation)
    this->d_cell_list.setLevelUp(0, 0, 0);

    // Assign the cell for Level Up (evaluate each particle)
    for (size_t i = 0; i < parEval.num; i++)
    {   
        // ... LEVEL 1 CHECK ...
        // Adaptation evaluation still based on vorticity value of the evaluated particle
        if (std::abs(parEval.vorticity[i]) >= 1.0e-5 * vor_max)  // Threshold set to be 1.0e-5 of the max value, \
                                                                    set lower than the particle redistribution
        {
            // ... LEVEL 2 CHECK ...
            // Do adaptation if the particle level != maxlevel (or lower) and the
            if (parEval.level[i] < Pars::max_level){
                // Assign the cell for Level Up
                this->d_cell_list.setLevelUp(parEval.basis_label[i], parEval.cell_label[i], 1);

                // Do the adaptation
                _adaptation = true;
            }
        }
    }
    
    // Adaptation procedure:
    // -> List all the cell to be devided (save the basis-cell-ID pair)
    // -> Find the neighboring cell to be divided (save the basis-cell-ID pair)
    // -> At each cell: divide the particle into the target level (only if par level < cell level)
    // -> Put the particle into the corresponding cell
    // -> Devide the cell
    
    // PROCEDURE 2 ! : Performing particle adaptation
    // *************
    if (_adaptation == true)
    {
        // Perform the adaptive algorithm
        this->d_cell_list.performAdaptive(baseParticle, 0, Pars::max_level, PARsize);

        // Particle adaptation summary time display
        t = clock() - t;
        printf("<+> Particle adaptation done ...\n");
        printf("<+> Number of particle after adaptation : %8d\n", baseParticle.num);
        printf("<-> Particle adaptation calculation \n");
        printf("    comp. time:                        [%f s]\n\n", (double)t/CLOCKS_PER_SEC);
        
        // Evaluate the verlet list of new particle distribution
        // this->neigbor_search(baseParticle);		// Already taken on the particle adaptation
    }else{
        printf("<+> No particle adaptation\n");
    }

    return _adaptation;
}

// ************************************************************************************
// ====================================================================================

// Particle redistribution
std::vector<std::vector<int>> neighbor::par2grid_neigbor_search(Particle& parEval, const Particle& parBase){
    // [1] Create the new list, Update the basis and cell ID of parEval
    // [2] Perform the neighbor search
    // Particle ID source: parEval
    // Particle ID neighbor target: parBase
    std::vector<std::vector<int>> ngh;
    return ngh;
}

