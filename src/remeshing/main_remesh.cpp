#include "remeshing.hpp"

// Particle Redistribution Main Manager
bool remeshing::get_remeshing(Particle &p, const Body &b, const int nIter)
{
	/*
	Inside Method:
	> Perform Distribution Adaptation
	  * Procedure:
	    - <Check whether the Adaptation is needed in the current state>
		- Evaluate the adapted particle (create a list)
		- Perform the split (available) and merging (still not available)
		- Update the cell List
		- Update the neighbor list
		- Update the particle base
	  * Needs:
	    - Determine the update region (Evaluate each particle vorticity value)
		- Needs:
		- Needs:
	
	> Perform Particle Redistribution (Must work no matter the Adaptation is not used)
	  * Procedure:
	    - Determine a new neighbor ID of the "Pattern Particle Base" - This is bitter when the adaptive performed
		- Perform the LSMPS B "Interpolate the vorticity of <Particle> for the <Pattern Particle Base>"
		- Change the <Particle> using <Pattern Particle Base>
		- Resize the other property (velocity, chi, etc.)
	  * Needs:
	    - Determine the interpolation neighbor (Take the neighbor of each neighbor particle ?)
		- Update the volume ratio in LSMPS B
		- 
	*/
	
	// computational time accumulation
    clock_t t;

	// ============================================================= //
	// ========== Particle Distribution Adaptation Update ========== //
	// ============================================================= //
    // Updating the temporary internal variable
	this->adap_flag = false;
	this->initial_num = this->particleBase.num;
	this->adtParSize = this->particleBase.s;
	
	// Threshold for adaptation -> Check adaptation only for 5 times
	if (Pars::flag_adaptive_dist && (nIter+1) % 5 == 0){
		this->adap_flag = true;
	}
	
	// Evaluating and Performing particle adaptation
	if (this->adap_flag == true){
		printf("\nPerform Particle Adaptation ...\n");
		// Perform the particle adaptation to the particle base
		printf("<+> Evaluation Particle Adaptation\n");
		
		// The current method will evaluate and peforming particle adaptation
		this->adap_flag = this->d_neighbor.particle_adaptation(p, this->particleBase, this->adtParSize);
	}

	// ============================================= //
	// ========== Particle Redistribution ========== //
	// ============================================= //
	t = clock();
	printf("\nPerform Particle Redistribution ...\n");
	// The package of particle redistribution
	this->redistribute_particles(p);

	// Update the size of the particle (since the particle adaptation still using the old particle size)
	// *** the code goes here ***
	if (this->adap_flag == true){
		this->d_neighbor.neigbor_search(p);
		this->particleBase.neighbor = p.neighbor;
	}
	
	// Particle redistribution summary time display
	t = clock() - t;
	printf("<+> Number of particle after remeshed   : %8d\n", p.num);
	printf("<-> Particle redistribution calculation \n");
	printf("    comp. time:                        [%f s]\n", (double)t/CLOCKS_PER_SEC);
	
	return this->adap_flag;
}

// Particle Redistribution Initiaization
void remeshing::remeshing_init(Particle &p){
	/* 
	The remeshing initialization procedure
	> Initialize the cell List
	> Neighbor search
	> Assign the particle base
	  * Particle Base Contains: pos, size, neighbor, cell ID, active
	*/
	printf("\n+------------ Particle Redistribution ------------+\n");

	// Generate the Cell List
    printf("Particle Redistribution ...\n");
    clock_t t = clock();

	// Initialize the cell List only for neighbor type 2
	if (Pars::opt_neighbor == 2){this->d_neighbor.cell_list_init(p);}

	// Evaluate the verlet list
	this->d_neighbor.neigbor_search(p);

	t = clock() - t;
	printf("\n\nTotal neighbor evaluation comp. time:  [%f s]\n\n", (double)t/CLOCKS_PER_SEC);
	
	// Assign the particle base
	this->particleBase = p;
	
	printf("+-------------------------------------------------+\n");
}