#include "geometry.hpp"

// Body Generation Manager
void geometry::generateBody(Body &b){
    printf("\n+------------------ Body Region ------------------+\n");

    // Computational time accumulation
    clock_t t = clock();

    // Create the solid body data according to body option in Pars namespace
    printf("Generating body surface node ...\n");
    if(Pars::DIM == 2 )
	{
		// Generate 2 Dimensional Body Geometry
		if (Pars::opt_body == 1) // Circular Cylinder [Subroutine CLEAR]
        {
			printf("<+> Type 1: Circular cylinder ... \n");
            this->cylinder_generator(b);
        }
		else if (Pars::opt_body == 2) // Flat Plate [Subroutine CLEAR]
        {
            printf("<+> Type 2: Flat plate ... \n");
		 	this->flat_plate_generator(b);
        }
		else if (Pars::opt_body == 3) // Plate (normal to freestream) [Subroutine CLEAR] 
        {
			printf("<+> Type 3: Normal plate ... \n");
            this->normal_plate_generator(b);
        }
		else if (Pars::opt_body == 4) // Squre Cylinder [Subroutine CLEAR] 
        {
			printf("<+> Type 4: Square cylinder ... \n");
            this->square_generator(b);
        }
        else if (Pars::opt_body == 5) // NACA Airfoil 
        {
            printf("<+> Type 5: NACA airfoil ... \n");
            this->naca_generator(b);
        }
		
	}
    else if(Pars::DIM == 3)
	{
		// Generate 3 Dimensional Body Geometry
        if (Pars::opt_body == 1) // Sphere 
        {
			printf("<+> Type 1: Sphere ... \n");
            this->ball_sphere(b);
        }
	}
    
    // Geometry computational time display
    t = clock() - t;
	printf("<-> solid body generation comp. time:  [%f s]\n", (double)t/CLOCKS_PER_SEC);

    // Geometry translational velocity
    // if (Pars::opt_u != 1){
        this->u_var(b);
    // }
    printf("+-------------------------------------------------+\n");
}

// Minimum Distance Calculation for one point
double geometry::distance_calc(std::vector<double> & pos, const Body& b){
	// The calculated variable
	double minDist;	// The minimum distance of particle to panel middle node
	double _dist;	// The temporarty distance calculation
	int panelNode;	// The panel node correspond to minimum distance
	double result;
	
	if (Pars::DIM == 2) // Calculation for 2D Space 
	{
		minDist = std::sqrt(std::pow(pos[0] - b.x_m[0],2) + 
				std::pow(pos[1] - b.y_m[0],2));
		panelNode = 0;
		for (int j = 1; j < b.n_panel; j++){
			_dist = std::sqrt(std::pow(pos[0] - b.x_m[j],2) + 
				std::pow(pos[1] - b.y_m[j],2));
			panelNode = minDist < _dist ? panelNode : j;
			minDist = minDist < _dist ? minDist : _dist;
		}
		
		// Inner product between minDist and normal direction
		result = (pos[0] - b.x_m[panelNode])*b.x_n[panelNode] + 
			     (pos[1] - b.y_m[panelNode])*b.y_n[panelNode];
	}
	else if (Pars::DIM == 3) // Calculation for 3D Space 
	{
		minDist = std::sqrt(std::pow(pos[0] - b.x_m[0],2) + 
				std::pow(pos[1] - b.y_m[0],2) + 
				std::pow(pos[2] - b.z_m[0],2));
		panelNode = 0;
		for (int j = 1; j < b.n_panel; j++){
			_dist = std::sqrt(std::pow(pos[0] - b.x_m[j],2) + 
				std::pow(pos[1] - b.y_m[j],2) + 
				std::pow(pos[2] - b.z_m[j],2));
			panelNode = minDist < _dist ? panelNode : j;
			minDist = minDist < _dist ? minDist : _dist;
		}
		
		// Inner product between minDist and normal direction
		result = (pos[0] - b.x_m[panelNode])*b.x_n[panelNode] + 
			 	 (pos[1] - b.y_m[panelNode])*b.y_n[panelNode] + 
				 (pos[2] - b.z_m[panelNode])*b.z_n[panelNode];
	}
	return(result);
}

// Minimum Distance Calculation for all particle
void geometry::distance_calc(Particle &p, const Body &b){
	// The calculated variable
	double minDist;	// The minimum distance of particle to panel middle node
	double _dist;	// The temporarty distance calculation
	int panelNode;	// The panel node correspond to minimum distance
	p.R.clear();	// Clear the particle distance to solid surface first
	p.R.resize(p.num,0.0);	// Update the size of distance
	if (Pars::DIM == 2) // Calculation for 2D Space 
	{
		for (int i = 0; i < p.num; i++){	// Iterate through all the particle inside the domain
			minDist = std::sqrt(std::pow(p.x[i] - b.x_m[0],2) + 
					std::pow(p.y[i] - b.y_m[0],2));
			panelNode = 0;
			for (int j = 1; j < b.n_panel; j++){
				_dist = std::sqrt(std::pow(p.x[i] - b.x_m[j],2) + 
					std::pow(p.y[i] - b.y_m[j],2));
				panelNode = minDist < _dist ? panelNode : j;
				minDist = minDist < _dist ? minDist : _dist;
			}
			
			// Inner product between minDist and normal direction
			p.R[i] = (p.x[i] - b.x_m[panelNode])*b.x_n[panelNode] + 
					  (p.y[i] - b.y_m[panelNode])*b.y_n[panelNode];
		}
	}
	else if (Pars::DIM == 3) // Calculation for 3D Space 
	{
		for (int i = 0; i < p.num; i++){	// Iterate through all the particle inside the domain
			minDist = std::sqrt(std::pow(p.x[i] - b.x_m[0],2) + 
					std::pow(p.y[i] - b.y_m[0],2) + 
					std::pow(p.z[i] - b.z_m[0],2));
			panelNode = 0;
			for (int j = 1; j < b.n_panel; j++){
				_dist = std::sqrt(std::pow(p.x[i] - b.x_m[j],2) + 
					std::pow(p.y[i] - b.y_m[j],2) + 
					std::pow(p.z[i] - b.z_m[j],2));
				panelNode = minDist < _dist ? panelNode : j;
				minDist = minDist < _dist ? minDist : _dist;
			}
			
			// Inner product between minDist and normal direction
			p.R[i] = ((p.x[i] - b.x_m[panelNode])*b.x_n[panelNode] + 
					  (p.y[i] - b.y_m[panelNode])*b.y_n[panelNode] + 
					  (p.z[i] - b.z_m[panelNode])*b.z_n[panelNode]);
		}
	}
}

/*
Note for geometry
*/