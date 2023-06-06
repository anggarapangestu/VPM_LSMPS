#ifndef INCLUDED_UTILS
#define INCLUDED_UTILS

#ifndef INCLUDED_GLOBAL
#include "global.hpp"
#endif

class Body
{
public:
	// NODE DATA
	int num;	// Number of Nodes
	// body coordinate
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

	// body extremes coordinate
	std::vector<double> min_pos;
	std::vector<double> max_pos;
	
	// body translational velocity
	std::vector<double> uT;
	std::vector<double> vT;
	std::vector<double> wT;
	
	// body rotational velocity
	std::vector<double> uR;
	std::vector<double> vR;
	std::vector<double> wR;
	
	// body deformation velocity
	std::vector<double> uDEF;
	std::vector<double> vDEF;
	std::vector<double> wDEF;
	
	// Body angular velocity (Still in 2D)
	std::vector<double> avelocity;

	// PANEL DATA
	int n_panel;	// Number of Panel
	// Panel middle coordinate (position)
	std::vector<double> x_m;
	std::vector<double> y_m;
	std::vector<double> z_m;

	// Panel middle coordinate
	std::vector<std::vector<int>> node_list;	// used only in 3D panel creation

	// Panel normal direction (directed inward the geometry) (vector)
	std::vector<double> x_n;
	std::vector<double> y_n;
	std::vector<double> z_n;
};

class Particle
{
public:
	// Basic Particle data
	int num;						// number of particle in the class
	std::vector<double> s;			// particle size (diameter)
	std::vector<std::vector<int>> neighbor; // list of neighbor index (verlet list add buffer gap length)

	// Data for Adaptive Tree Cell
	std::vector<int> basis_label;	// The label of basis cell ID
	std::vector<int> cell_label;	// The label of cell ID
	std::vector<int> level;			// The finenest level of particle resolution
	
	// Data for Marking active particle
	std::vector<bool> isActive;		// indicates active particles [@param: transient particle]
	std::vector<int> isActiveIndex; // indicates active particles' index [@param: transient particle] <???>

	// particle coordinates (position)
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	
 	// particle velocities and vorticities data
	std::vector<double> u;
	std::vector<double> v;
	std::vector<double> w;
	std::vector<double> gz;				// vortex strength (gamma)
	std::vector<double> vorticity;		// The scalar field
	std::vector<double> P;				// The scalar field of pressure
	
	// Penalization parameters
	std::vector<double> chi;
	std::vector<double> R;		// The particle shortest distance from the solid surface (negative inside, positive outside)
	
	// Additional
	std::vector<int> hash_cell;	// STILL CHECKED

	// Boundary terms
	std::vector<int> isboundary;		// <???> Maybe good for boundary treatment
	std::vector<bool> inside; 
	std::vector<double> boundaryval;
};

// Tambahan untuk hashing 3-d: Ical
// Need to be deleted
class Cell
{
public:
    std::vector<int> x; 
    std::vector<int> y;
    std::vector<int> z;
    std::vector<int> num_grid; 
    std::vector<std::vector<int>> particle_inside; 
};

#endif