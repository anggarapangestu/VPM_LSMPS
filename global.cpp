#include "global.hpp"

namespace Pars
{
// Predescription for each parameter variable
// [A] -> Already given in this program
// [I] -> Input for PARAMETER form user
// [O] -> Input for OPTION form user, if not will use default
// [C] -> Calculated form user input
    
    // ================== Parameter for Testing ================== //
    extern const int testtype = 1;             // 1 = impulsively, 2 = flowing

    // ================== Basic Math Parameter/Constant ================== //
    extern const double pi = M_PI /*3.14159265358979323846*/;
    extern const bool flag_save_body = true;       // The flag to save body data
    extern const bool flag_save_cell = false;       // The flag to save cell list data
    extern const bool flag_save_parameter = true;  // The flag to save simulation parameter
    extern const bool flag_save_sim_time = true;   // The flag to save simulation time and computational time
    extern const bool flag_disp_stability = true;  // The flag to display stability parameter at each iteration
    extern const bool flag_pressure_calc = true;   // The flag to calculate pressure
    extern const bool flag_adaptive_dist = true;   // The flag for adaptive particle distribution

    extern const int opt_cont = 0;                 // Option for starting simulation 0:= Start at initial, 1:= continue the stoped running

    // ================== Simulation Fundamental Option ================== //
    extern const int DIM = 2;       // [I] Spatial domain dimension; 2= 2D simulation, 3= 3D simulation
    extern const int iterative = 1;     // 1 = classical brinkmann, 2 = iterative brinkmann
    extern const double lxdom = 6.0e0;         // Initial condition domain length of particle distribution on x-axis
    extern const double lydom = 5.0e0;         // Initial condition domain length of particle distribution on y-axis
    extern const double lzdom = 0.0e0;         // Initial domain length (z-axis)
    extern const double xdom = 2.0e0;          // negative x-direction domain length 
    extern const double xcenter = 0.0e0; // Initial domain center (x-axis), default 0.0
    extern const double ycenter = 0.0e0; // Initial domain center (y-axis), default 0.0

    /*domain illustration
    |------------------------------|
    |                              |
    |                              |
    |     (xcenter, ycenter)       |
    |<------------>*               | lydom
    |     xdom                     |
    |                              |
    |                              |
    |------------------------------|
                lxdom
    */

    // Initialization
    extern const int init_type = 1;            // 0= random ditribution, 1= uniform distribution
    extern const double dev = 0.2;             // [0.2, 0.25] particle deivation in irregular particle distribution
    extern const double c = 1.0;               // [0.5, 0.9, 1.4] overlapping ration

    // ================== Simulation Solid Parameter ================== //
    extern const int opt_body = 1;  // [O] The body type option; \
                                            1:= 2D - Circular cylinder, 3D - Sphere;\
                                            2:= Flat plate;\
                                            3:= Normal plate;\
                                            4:= Square cylinder;\
                                            5:= NACA airfoil;
    extern const int opt_init = 4;  // [O] The initialization type option; \
                                            0:= Testing Resolution;\
                                            1:= Single Resolution;\
                                            2:= Multi-Res Single Block;\
                                            3:= Multi-Res Multi Block; // [Not Ready YET] \
                                            4:= Multi-Res Body Adjusted;
    int opt_neighbor = 2;  // [O] The neighbor evaluation type option; \
                                            0:= Direct Neighbor Search;\
                                            1:= Linked List;\
                                            2:= Cell List;\
                                            3:= Spatial Hash;
    
    extern const int n_a = 400;       // [I] Number of body surface node
    // Body velocity
    extern const int opt_u = 1;     // [O] Body velocity option; 1= Suddenly start and stop (Instanteneosly appear), 2= smoothly start and stop, 3= moving
    extern const double ubody = 0.0e0;  // [I] Body x-direction velocity
    extern const double vbody = 0.0e0;  // [I] Body y-direction velocity
    // Whole Geometry Parameter
    // NACA4-digit series: NACA(m*100)(p*10)(t*100) 2412,4612, Cylinder, or Flat plate
    extern const double m_a = 0.0e0;    // [I] max thickness
    extern const double p_a = 0.0e0;    // [I] maximum thickness position in percent of chord
    extern const double t_a = 0.15e0;   // [I] is maximum thickness of airfoil in tenth of chord
    extern const double acx = 0.25*1;   // [I] aerodynamic center location
    extern const double acy = 0.0e0;    // [I] aerodynamic center location
    extern const double alpha_a = 0.0e0;// [I] alpha in degree , rotated angle (attached angle), ! Rotated airfoil
    extern const double H_star = 1.0e0 / 10.0e0; // [I] H*=H/L
    extern const double ly = 2.0e0;     // [I] flat plate: length, sphere: diameter (flat plate: length, sphere: diameter, Square: side.)
    extern const double lx = 2.0e0;     // [I] sphere: diameter
    // chord*length ! lx*ly ! Flatplate, if sphere  A= pi*Df^2/4
    extern const double Df = 1.0e0;     // [I] Reference length, flat plate: chord, sphere: diameter

    // ================== Incompressible Fluid Properties ================== //
    extern const double RE = 550.0e0;  // [I] The reynold number (U)
    extern const double u_inf = 1.0e0;  // [I] Freestream x-direction velocity
    extern const double v_inf = 0.0e0;  // [I] Freestream y-direction velocity
    extern const double w_inf = 0.0e0;  // [I] Freestream z-direction velocity
    extern const double RHO = 1.0e0;    // [I] Fluid density
    extern const double U_inf =         // [C] Freestream velocity magnitude
        std::sqrt(u_inf*u_inf + v_inf*v_inf);
    extern const double NU = U_inf * Df / RE;   // [C] Kinematic viscousity
    extern const double MU = RHO * NU;          // [C] Dynamic viscousity

    // ================== Computation Parameter ================== //
    extern const double sigma = 0.01e0;    // [I] Particle core size; ? default: 0.0025, dt~(1/2)*sigma - (1/3)*sigma
    extern const double dt = 0.01e0;       // [I] Simulation time step; ? default:0.001, dt <= phi_s*sigma^2/vis (Ploumhans [2000]) OR dt = phi_s*sigma^2*Courant^2/vis, where 0 < Courant <= 1
    extern const double phi_s = 0.595e0;    // [A] (Ploumhans [2000]) for the Euler explicit scheme,phi_s = 0.595, !for the Adams–Bashforth 2 scheme, phi_s = 0.297, BUt can not use AB2 due to redistribution changes np--> need Biot-Savart again
    extern const double Courant = U_inf * dt / sigma;           // [C] Courant number(C): C = U_inf * dt / sigma, where 0 < C <= 1
    extern const double Diffusion = NU * dt / (sigma*sigma);    // [C] Diffusion number(Phi): Phi =  NU * dt / sigma^2, where 0 < Phi <= 0.5
    // extern const double Courant = std::sqrt(dt * vis / phi_s) / sigma;   // sigma = sqrt(dt*vis/phi_s)/Courant, where 0 < Courant <= 1
    extern const double r_sup = 3.2;        // [I] Neighbor evaluation support domain ratio size
    extern const double r_buff = 1.2;       // [I] The buffer region radius factor
    extern const int max_level = 2;         // [I] Number of resolution level
                                                // -> Level is counted from 0 as the larger particle size
                                                // -> Maximum level is the finest particle size
                                                // -> The number of resolution is 'max_level' + 1

    // ================== Penalization Parameter ================== //
    extern const int opt_pen = 2;             // [O] The penalization type option; \
                                                    1:= implicit penalization; (lambda = 1.0e4) \
                                                    2:= semi-implicit penalization,Rasmussen 2011 mixed with explicit; \
                                                    3:= Explicit scheme, Rasmussen 2011 (lambda = 1.0d0/dt)
    int opt_kaipen = 0;                       // [O] The χ (chi/kai) option; 0:= no recalculate kai; 1:= recalculated by kai~ ! to avoid singularities Rasmussen phD 2011
    extern const int numpen = 6;             // [I] Number of points from surface for penalization domain;beberapa referensi ngasihnya 6
    // Variable r_epsilon
    extern const double hmollif = 1.0e0 * std::sqrt(2.0e0) * sigma;// / 2.0e0;    // [I] The mollification length [2sqrt(2)*sigma] <- 2D, [sqrt(2)*sigma]
    // // Constant r_epsilon
    // extern const double hmollif = 0.01 * Df;// / 2.0e0;    // [I] The mollification length [2sqrt(2)*sigma] <- 2D, [sqrt(2)*sigma]
    double lambda = 1.0e4;                        // [C] Value of penalization

    extern const bool opt_run_simulation = true; // true = save data during penalization; vice versa.

    // ================== Saving Parameter ================== //
    extern const double simulation_time = 10.0e0;    // The simulation time??
    extern const double comtime_sf = 20.0e0 * 500.0e0 * 600.0e0;   // % Save file frequently after how many [second]=[hour]*60*60, only in case of running is stopped suddenly, but nt_sf take long days
    extern const int force_type = 1;        // 0 no force saved; 1 for penalization; 2 for NOCA
    extern const int opt_extra_data = 3;                           // = 1 not print extra data, print only common data
                                                                // = 2 print common data and extra data at the same time,
                                                                // = 3 print extra data from common data, MUST run case 1/2 first
    // in case of opt_extra_data = 1/2:
    extern const int nt = std::ceil(simulation_time / dt); // number of iterations
    extern const int nt_data = 10;                    // number of data saved, added 16:04 March 23 by Angga
    extern const int step_inv = 20;//std::ceil(nt / nt_data);   // step interval for data saving
    extern const int nt_sf = 50;                    // save data per ( nt_sf ) time step, for saving storage MEMORY and for case of running is stopped suddenly

    // ================== Other New Parameter ================== //
    extern const double init_act_dist = 0.3;           // The initial active distance factor toward particle size
    // extern const double simulation_time;    // 


extern const double r_sigma = std::sqrt(sigma/pi);
extern const double alpha = 0.8e0;                                   // (Overlaping factor)
extern const double gamtrim = 1.0e-04;                               // parameter for Extruding blobs, ! 2D: 1.0e-04 Ploumhans 2000, 3D 1.0e-04 (Ploumhans [2002])
extern const double re_trsh = 1.0e-04;                               // Re_h,trsh = 1.0e-04 (Ploumhans [2002]) 3D
extern const double area  = pi*std::pow(Df,2)/4;                     // chord*length ! lz*ly ! Flatplate, if sphere  A= pi*Df^2/4

extern const double x_new_1 = -5; 
extern const double x_new_2 = 30;
extern const double y_new_1 = -10;
extern const double y_new_2 = 10;
//Note: Final shape of multiblock domain will be following this formula  side x = ( lxdom + (2*10*Pars:sigma)) - 2*Pars::sigma. Y side will also adopt the formula


// Multiresolution parameters
extern const double D_0 = 2 * sigma;              // maximum core size:::4
extern const double dcmin = 0.4;                  // minimum equivalent distance

//Spatial Hashing Parameters 
extern const double xmin = -2; 
extern const double ymin = -2;
extern const double zmin = -2;
extern const double xmax = 2; 
extern const double ymax = 2;
extern const double zmax = 2;
extern const double cell_size = 0.25 * (xmax - xmin); // !!!!size of cell!!!!!


// DC PSE operators parameter ===================================== //
extern const double beta = 2;
extern const double epsilon = c * sigma;
extern const int l_2_0_2 = 9;                 // moment condition for Q(2,0) or Q(0,2) with 2nd order accuracy
extern const int l_2_0_4 = 20;                // moment condition for Q(2,0) or Q(0,2) with 4th order accuracy
extern const int l_1_0_2 = 6;                 // moment condition for Q(1,0) or Q(0,1) with 2nd order accuracy
extern const int l_1_0_4 = 15;                // moment condition for Q(1,0) or Q(0,1) with 4th order accuracy
extern const int l_0_0_3 = 6;                 // moment condition for Q(0,0) with 3rd order accuracy

// in case of continue

extern const int it_stop_full = 60166;         // the last file which have full data (before/at stop) (last file saved by comtime_sf or nt_sf), [changeable]

// in case of opt_extra_data = 3 print extra data from common data, MUST run case 1/2 first:
extern const int it_start = 0;                 // the first step you want to get extra data for postprocessing
extern const int it_end = 101 /*nt*/;          // the last step you want to get extra data for postprocessing
// in case of opt_extra_data = 2/3:
extern const int opt_extra_pres = 0;           // = 0/1  not/ print pressure on grid
extern const int opt_extra_pres_wall = 0;      // = 0/1  not/ print pressure on surface of body

// ====== BOX Generation for post-processing extra data =========== //
extern const int opt_gdata_vel = 2;            // =1 Using Biot-savart calculate Velocity on grid via Vorticity Strength on grid, =2 Remeshing Velocity from particle velocity

// parameter for link list ======================================== //
extern const int dim = 2;                      // dimension
extern const int skf = 2;                      // Gauss spline
extern const int r_scale = 3;                  // For box size. 

// Tree option
extern const double expTree = 1.0;  // The percentage of tree root cell size expansion

// FMM option
extern const int P_max = 15;          // Number of the expansion
extern const int level_max = 7;       // The tree level of the finest cell
extern const int level_min = 6;       // The tree level of the coarse cell
extern const int n_max = 40;          // Maximum number of the particle inside the cell
extern const double tol_dist = 1.0e-10;   // The maximum tolerance of cell size

// parameter for FMM ============================================== //
extern const int lmax = 7;            // max of FMM levels
extern const int npcm = 100;          // Number of pole expansion
extern const int nfwrd = 16;
// std::pow(4,(lmax+1));   //nfwrd * (std::pow(4, (lmax + 1)) - 4) / 3;
extern const int nbmax = 4 * (std::pow(4,(lmax))-1) / (4-1);   // Geometric Series
extern const int nbl1 = 8 * lmax;
extern const int nbl2 = 256;
extern const int nbl3 = 256;
extern const int nbl4 = 16;
extern const int nbmrl = nbmax + 5; //nfwrd * 4096; //!!need to change if so that nbmrl >= 4^lmax
extern const double tol2 = 1.0e-15;

// FMM and FMMF option and their box size for solution ============ //
extern const int ioptfmm = 1;                   // = 0 Direct, = 1 FMM
extern const int n_s = 50;                      // 10 if using fmmo *******
extern const int ndp = 10;                      // number of expansion
extern const int icutoff = 2;                   // 0.singular ;  1.super (high-oder) algebraic  ; 2. Gaussian  ; 3 super Gaussian
                                                // if par_diff = 0 (PSE): take 2or3; if ioptfmm =1(FMM): take 1; if ioptfmm =1(FMMF): take 2or3
extern const int iopt_inter = 1;                //
extern const int par_ext = 1;                   // =1(external flow case), =0(internal flow case)

// REMESHING, DIFFUSION, SEARCHING, CONVECTION, STRETCHING, LES === //
extern const int nrmsh = 1;//1;                 // Remeshing every "nrmsh" timestep
extern const int opt_remesh_W = 1;              // =1 using M4', =2 Using M6, =3 M4 isotropic (Chatelain,P [2005]), =4 P14
extern const int par_split = 5;                 // 4,5,7,9 only, otherwise just spreading core size
extern const int parsplit = 4;                  //
extern const int xlimit = 10000;                // in order to eliminate far field splitting
extern const int opt_sptadapt = 1;              // = 1 using spatial_adaptation for CSM, =0 withou spatial_adaptation
extern const int opt_search = 2;                // =1 using direct searching O(N^2); =2 using link_list O(N)
extern const int it_start_les = 2;              // if  it_start_les =1, with zero gpx gpy gpz then velocity gradient dudx dudy ... = 0 at it =1, then strain(it =1)=0...=> Cr2 = NaN
extern const int diffusion_opt = 0;             // PSE scheme. 0=original PSE, 1=DC-PSE

// VIBRATION / STRUCTURAL PROPERTIES : !!STILL IN PROGRESS, RESULTS NOT GOOD ENOUGH!!
extern const int vib = 0;                       // =0 no vibration simulation, =1 circular cylinder simulation, =2 flow induced pitch oscillations
double SpringConst = 0.0;                // Spring constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
double DamperConst = 0.0;                // Damper constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
double mass = 0.0;                       //Object mass :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
double inertia = 0.0;                    // Inertia :: !!STILL UNUSED!!
//extern double m_d = dens*area;                // !! rho times volume enclosed !! 
double m_d = 0.0;                        // !! rho times volume enclosed !! UNUSED, SET TO 0.0 
double tetha = 0.0;                      //!! ROTATION PARAMETER, STILL BUGGED
double tetha_dot = 0.0;                  //!! ROTATION PARAMETER, STILL BUGGED
extern const double m_star = 2.5;               // NONDIMENSIONAL MASS
extern const double k_star = 4.96;              // NONDIMENSIONAL SPRING CONST
extern const double c_star = 0;                 //NONDIMENSIONAL DAMPER CONST
extern const double t_star;                     //NONDIMENSIONAL TIME VARIABLE
double gaya = 0.0;                       // !! FORCE INITIALIZATION !!
double momen = 0.0;                      // !! MOMENT INITIALIZATION !!
extern const double chi_star = 0.15;            //!! ROTATION PARAMETER, STILL BUGGED
extern const double i_star = 4.1;               //!! ROTATION PARAMETER, STILL BUGGED
extern const double U_star = 6.6;               //!! ROTATION PARAMETER, STILL BUGGED
extern const double tetha_nol = 15.0e0;         //equilibrium angular position of spring [deg] //!! ROTATION PARAMETER, STILL BUGGED


namespace geom
{
extern const int edge = 2;
}

} // namespace Pars