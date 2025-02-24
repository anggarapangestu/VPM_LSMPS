#ifndef INCLUDED_UTILS
#define INCLUDED_UTILS

#include "global.hpp"

/**
 *  @brief A property set of single body made of node or panel.
 *
 *  @headerfile Utils.hpp
 */
struct Body{
	// GENERAL DATA
	// ************
	// Body extremes coordinate
	double cen_pos[DIM];		// Body center coordinate
	double min_pos[DIM];		// Minimum extreme in each basis dimension
	double max_pos[DIM];		// Maximum extreme in each basis dimension

	// Body node translational velocity (*at each simulation iteration)
	std::vector<double> uT;		// Translation velocity in x direction
	std::vector<double> vT;		// Translation velocity in y direction
	std::vector<double> wT;		// Translation velocity in z direction
	
	// // body rotational velocity (*at each simulation iteration <?>)
	// std::vector<double> uR;
	// std::vector<double> vR;
	// std::vector<double> wR;
	
	// // body deformation velocity (*at each body note point <?>)
	// std::vector<double> uDEF;
	// std::vector<double> vDEF;
	// std::vector<double> wDEF;
	
	// // Body angular velocity (Still in 2D)
	// std::vector<double> avelocity;

	// NODE DATA
	// *********
	int n_node;					// Number of Nodes

	// Body node coordinate
	std::vector<double> x;		// Node coordinate in x basis
	std::vector<double> y;		// Node coordinate in y basis
	std::vector<double> z;		// Node coordinate in z basis

	// PANEL DATA
	// **********
	int n_panel;				// Number of Panel
	
	// Body panel middle coordinate (position)
	std::vector<double> x_m;	// Panel Coordinate in x basis
	std::vector<double> y_m;	// Coordinate in x basis
	std::vector<double> z_m;	// Coordinate in x basis
	std::vector<double> size;	// The panel size (2D: length, 3D: area)
	
	// IMPORTANT: normal vector directed outward the body (pointing outside)
	// Panel normal vector
	std::vector<double> x_n;
	std::vector<double> y_n;
	std::vector<double> z_n;

	// // Panel middle coordinate
	// std::vector<std::vector<int>> node_list;	// used only in 3D panel creation

	// Constructor
	Body(){};

	// Deconstructor
	~Body(){};
};

/**
 *  @brief A data set of particle inside the domain.
 *
 *  @headerfile Utils.hpp
 */
struct Particle{
	// Physical Properties
	// *******************
	// Particle position coordinates
	std::vector<double> x;			// Coordinate in x basis
	std::vector<double> y;			// Coordinate in y basis
	std::vector<double> z;			// Coordinate in z basis
	
 	// Particle physical properties
	std::vector<double> u;			// Velocity in x direction
	std::vector<double> v;			// Velocity in y direction
	std::vector<double> w;			// Velocity in z direction
	std::vector<double> gx;			// Vortex strength x (gamma_x) [*]
	std::vector<double> gy;			// Vortex strength y (gamma_y) [*]
	std::vector<double> gz;			// Vortex strength z (gamma_z) [*]
	std::vector<double> vortx;      // Vorticity x (omega_x)
	std::vector<double> vorty;      // Vorticity y (omega_y)
	std::vector<double> vortz;      // Vorticity z (omega_z)
	std::vector<double> vorticity;	// Vorticity scalar field (or Vorticity absolute value)
	std::vector<double> P;			// Pressure scalar field
	// Note on [*] : Not gonna used, calculate once and use once in simulation, seem redundant on memory
	
	// Interpolation Additional
	// ************************
	// For data interpolation
	std::vector<double> Q;      	// The Q criterion value : The second invariant of velocity gradient
	std::vector<double> L2;      	// The L2 criterion value : median of (S^2 + W^2) tensor eigenvalue (Jeong,1995)

	// Testing Properties
	// ******************
	// For perlman vorticity neumann
	std::vector<double> dudx;       // Velocity x differential toward x
	std::vector<double> dudy;       // Velocity x differential toward y
	std::vector<double> dvdx;       // Velocity y differential toward x
	std::vector<double> dvdy;       // Velocity y differential toward y

	// A testing data for poisson calculation
	std::vector<double> phi_n;      // The numerical poisson solution 
	std::vector<double> phi_a;		// The analytical poisson solution
	std::vector<double> F;			// The poisson source

	// A testing data for vorticity calculation
	std::vector<double> vort_a;		// Analytical vorticity
	std::vector<double> u_a;		// Analytical velocity in x direction
	std::vector<double> v_a;		// Analytical velocity in y direction
	std::vector<double> w_a;		// Analytical velocity in z direction

	// A bypass value for vorticity differentials (For Adaptation Error prediction)
	std::vector<double> absLapVortX;	// The absolute laplacian of X vorticicy
	std::vector<double> absLapVortY;	// The absolute laplacian of Y vorticicy
	std::vector<double> absLapVortZ;	// The absolute laplacian of Z vorticicy


	// For a test data
	// ... Please add the test data

	// Computational Properties
	// ************************
	// Basic computational Properties
	int num;						// Number of particle in the class
	int currStep;					// The current iteration step
	bool isAdapt;					// A flag to said particle is adapted
	std::vector<double> s;			// Particle size (diameter)
	std::vector<int> level;			// The level of particle resolution (related to size)
	std::vector<std::vector<int>> neighbor; // Index list of neighbor particle [A verlet list]

	// [GROUP 1] Data using *Grid Node* container
	std::vector<int> nodeID;		// ID label of node container

	// [GROUP 2] Data using *Adaptive Tree Cell* container
	std::vector<int> basis_label;	// ID label of basis cell
	std::vector<int> cell_label;	// ID label of cell

	// Penalization parameters
	std::vector<double> chi;		// Particle penalization mask value
	std::vector<double> R;			// Shortest distance from the solid surface; sign: [-] inward, [+] outward

	// Flag Marker and near body variable
	std::vector<int> bodyPart;			// Nearest body part ID [-1:= for not near to any body part]
	std::vector<bool> isActive;			// Indicates active particle => Particle that have a vorticity value [@param: transient particle]
	std::vector<bool> isNearSurface;	// Flag for near body surface
	std::vector<bool> insideBody;		// Flag for inside body

	/* Near Body Criteria Illustration
       	 __________________________________________________
		|        .             x _____x_______   x     .   |       Note:
		| .          .  *    x  |#2 x   x     \       .    |		> In the given illustration 3 objects is existing
		|     *_________        |      x       \ x   .     |		   inside the domain
		|   * /#1     * \ *   x |____________x__\  x     . |		> Particle is grouped in 4 criteria:
		|    /  *  *     \        x           x            |		   - [.] : free particle (not near with any body)
		|   /  *     *   / *   x  o  ____o__________o__ o  |           - [*] : particle near and inside body #1
		|   \     *   * /       o   |#3 o      o o     |   |           - [x] : particle near and inside body #2
		|  * \_________/ *        o |               o  |   |           - [o] : particle near and inside body #3
		|   *                    o  | o    o o    o    |   |        > The "bodyPart" store the nearest body part of 
		|     *         .      .    |   o      o     o | o |           the current particle, for free particle group
		|    .    .                o|__________________|o  |           value is set to -1
		|  .  .           .    .     o   o    o     o      |           
		|__________________________________________________|           
	*/

	// // Boundary terms <?> For boundary treatment <?> Not really need for LSMPS
	std::vector<int> boundaryLoc;		// The location of boundary [-2: left, -1: bottom, 0: inside, 1:up, 2: right]
	std::vector<double> boundaryVal;	// The dirichlet value of boundary condition
	// std::vector<bool> inside;

	// Data for LSMPS
	// std::vector<double> F;
	std::vector<double> Fx;
	std::vector<double> Fy;
	std::vector<double> Fx2;
	std::vector<double> Fy2;

	// Analitical Data
	std::vector<double> F_a;
	std::vector<double> Fx_a;
	std::vector<double> Fy_a;
	std::vector<double> Fx2_a;
	std::vector<double> Fy2_a;

	// Constructor
	Particle() : num(0), currStep(0), isAdapt(false){};

	// Deconstructor
	~Particle(){};
};

/**
 *  @brief A simulation utilities package list.
 *
 *  @headerfile Utils.hpp
 */
class simUtil{
private:
	int iterDigitLen;
	double E_0;		// The initial enstrophy
	double J_0;		// The initial angular momentum
	double G_0;		// The initial circulation
public:
	// Displaying and utilities calculation

	void startCounter(int step);
	void printHeader(int step);
	std::string saveName(const int step);
	void predictCompTime(int step, double _currTime);

	// Reduce particle size
	void trimParticleDomain(Particle &_particle);

	// Evaluation Function
	void saveResidual(Particle &_particle, int step);
	// void saveCompTime(int step, double _currTime);
	void addVorPertubation(Particle &_particle);

	// Physical properties calculation
	void calculate_save_enstrophy(Particle &_particle, int _step);
	void calculate_save_kinetic(Particle &_particle, int _step);
	void calculate_save_dissipation(Particle &par, int it);

	double calculate_L2_Norm(const Particle &_particle, const std::vector<double> &f_n, const std::vector<double> &f_a);
	void calculate_L_Norm(std::vector<double>& Norm, const Particle &_particle, const std::vector<double> &f_n, const std::vector<double> &f_a);
	void velocity_L2Norm(Particle &_particle, int _step, int _func, int _type);
	void vorticity_L2Norm(Particle &_particle, int _step, int _func);
	
	void timingCount(double duration, std::string nameTime, int end);
};

#endif