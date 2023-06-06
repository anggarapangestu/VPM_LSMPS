#ifndef INCLUDED_GLOBAL
#define INCLUDED_GLOBAL

// Mains and basics
#include <iostream>
#include <cmath>
#include <vector>
#include <time.h>   // clock_t
#include <chrono>   // clock_t

// Additionals
#include <math.h>
#include <omp.h>    // pragma omp ...
#include <complex>  // std::complex<>
#include <numeric>
#include <algorithm>

namespace Pars
{
// Predescription for each parameter variable
// [A] -> Already given in this program (default constant)
// [I] -> Input for PARAMETER form user 
// [O] -> Input for OPTION form user, if not will use default
// [C] -> Calculated form user input
    
    // ================== Parameter for Testing ================== //
    extern const int testtype;  

    // ================== Basic Math Parameter/Constant ================== //
    extern const double pi;     // [A] Math pi

    // ================== Simulation flag option ================== //
    extern const bool flag_save_body;       // The flag to save body data
    extern const bool flag_save_cell;       // The flag to save cell list data
    extern const bool flag_save_parameter;  // The flag to save simulation parameter
    extern const bool flag_save_sim_time;   // The flag to save simulation time and computational time
    extern const bool flag_disp_stability;  // The flag to display stability parameter at each iteration
    extern const bool flag_pressure_calc;   // The flag to calculate pressure
    extern const bool flag_adaptive_dist;   // The flag for adaptive particle distribution

    extern const int opt_cont;     // =0/1 Not/continue the stoped running

    // ================== Simulation Fundamental Option ================== //
    extern const int DIM;       // [I] Spatial domain dimension; 2= 2D simulation, 3= 3D simulation
    extern const double lxdom;      // Initial domain length (x-axis)
    extern const double lydom;      // Initial domain length (y-axis)
    extern const double lzdom;      // Initial domain length (z-axis)
    extern const double xdom;
    extern const double xcenter;    // Initial domain center (x-axis), default 0.0
    extern const double ycenter;    // Initial domain center (y-axis),  default 0.0

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
    extern const int init_type;
    extern const double dev;
    extern const double c;

    // ================== Incompressible Fluid Properties ================== //
    extern const double RE;     // [I] The reynold number (U)
    extern const double u_inf;  // [I] Freestream x-direction velocity
    extern const double v_inf;  // [I] Freestream y-direction velocity
    extern const double w_inf;  // [I] Freestream z-direction velocity
    extern const double RHO;    // [I] Fluid density
    extern const double U_inf;  // [C] Freestream velocity magnitude
    extern const double NU;     // [C] Kinematic viscousity
    extern const double MU;     // [C] Dynamic viscousity

    // ================== Computation Parameter ================== //
    extern const double sigma;      // [I] Particle core size; ? default: 0.0025, dt~(1/2)*sigma - (1/3)*sigma
    extern const double dt;         // [I] Simulation time step; ? default:0.001, dt <= phi_s*sigma^2/vis (Ploumhans [2000]) OR dt = phi_s*sigma^2*Courant^2/vis, where 0 < Courant <= 1
    extern const double phi_s;      // [A] (Ploumhans [2000]) for the Euler explicit scheme,phi_s = 0.595, !for the Adams–Bashforth 2 scheme, phi_s = 0.297, BUt can not use AB2 due to redistribution changes np--> need Biot-Savart again
    extern const double Courant;    // [C] Courant number(C): C = U_inf * dt / sigma, where 0 < C <= 1
    extern const double Diffusion;  // [C] Diffusion number(Phi): Phi =  NU * dt / sigma^2, where 0 < Phi <= 0.5
    extern const double r_sup;      // [I] Neighbor evaluation support domain ratio size
    extern const double r_buff;     // [I] The buffer region radius factor
    extern const int max_level;     // [I] Number of resolution level
                                        // -> Level is counted from 0 as the larger particle size
                                        // -> Maximum level is the finest particle size
                                        // -> The default max level is 0 (which is single resolution)
    
    // ================== Simulation Solid Parameter ================== //
    extern const double Df;     // [I] Reference length, flat plate: chord, sphere: diameter
    extern const int opt_body;  // [O] The body type option
    extern const int opt_init;  // [O] The initialization type option
    extern const int n_a;       // [I] Number of body surface node
    extern const double H_star; // [I] H*=W/L the width to length ratio
    // Body velocity
    extern const int opt_u;     // [O] Body velocity option; 1= Suddenly start and stop (Instanteneosly appear), 2= smoothly start and stop, 3= moving
    extern const double ubody;  // [I] Body x-direction velocity
    extern const double vbody;  // [I] Body y-direction velocity
    // Whole Geometry Parameter
    // NACA4-digit series: NACA(m*100)(p*10)(t*100) 2412,4612, Cylinder, or Flat plate
    extern const double m_a;    // [I] max thickness
    extern const double p_a;    // [I] maximum thickness position in percent of chord
    extern const double t_a;    // [I] is maximum thickness of airfoil in tenth of chord
    extern const double acx ;   // [I] aerodynamic center location
    extern const double acy ;   // [I] aerodynamic center location
    extern const double alpha_a;// [I] alpha in degree , rotated angle (attached angle), ! Rotated airfoil
    extern const double ly;     // [I] flat plate: length, sphere: diameter
    extern const double lx;     // [I] sphere: diameter
    
    // ================== Penalization Parameter ================== //
    extern const int opt_pen;       // [O] The penalization type option; 1:= semi-implicit penalization, Rasmussen 2011 mixed with explicit; 2:= Explicit scheme, Rasmussen 2011
    extern const int opt_lampen;    // [O] The lambda option; 1:= lambda = 1.0d4; 2:= lambda = 1.0d0/dt
    extern int opt_kaipen;          // [O] The χ (chi/kai) option; 0:= no recalculate kai; 1:= recalculated by kai~ ! to avoid singularities Rasmussen phD 2011
    extern const int numpen;        // [I] Number of offset particle from body surface for penalization domain
    extern const double hmollif;    // [I] The mollification length
    extern double lambda;           // [C] Value of penalization
    extern const int iterative;     // A2D

    extern const bool opt_run_simulation; // true = save data during penalization; vice versa.

    // ================== Saving Parameter ================== //
    extern const double simulation_time;    // The simulation time??
    extern const double comtime_sf;         // % Save file frequently after how many [second]=[hour]*60*60, only in case of running is stopped suddenly, but nt_sf take long days
    extern const int opt_extra_data;        // 
    extern const int force_type;            // asd
    extern const int nt;         // number of iterations step
    extern const int nt_data;    // number of data saved, added 16:04 March 23 by Angga
    extern const int step_inv;   // step interval for data saving
    extern const int nt_sf; // save data per ( nt_sf ) time step, for saving storage MEMORY and for case of running is stopped suddenly

    // ================== Other New Parameter ================== //
    extern const double init_act_dist;           // The initial active distance factor toward particle size
    // extern const double simulation_time;    // 

    // in case of opt_extra_data = 1/2:
    // in case of continue
    extern const int it_stop_full; // the last file which have full data (before/at stop) (last file saved by comtime_sf or nt_sf), [changeable]
    // in case of opt_extra_data = 3 print extra data from common data, MUST run case 1/2 first:
    extern const int it_start; // the first step you want to get extra data for postprocessing
    extern const int it_end;   // the last step you want to get extra data for postprocessing
    // in case of opt_extra_data = 2/3:
    extern const int opt_extra_pres;      // = 0/1  not/ print pressure on grid
    extern const int opt_extra_pres_wall; // = 0/1  not/ print pressure on surface of body

    // THE BELOW PARAMETER IS IN A CONSIDERATION !!!

    //================== Simulation Selection ==========================//
    extern const int sim;
    extern const double alpha;   // (Overlaping factor)
    extern const double gamtrim; // parameter for Extruding blobs, ! 2D: 1.0e-04 Ploumhans 2000, 3D 1.0e-04 (Ploumhans [2002])
    extern const double re_trsh; // Re_h,trsh = 1.0e-04 (Ploumhans [2002]) 3D
    extern const double area;    // chord*length ! lz*ly ! Flatplate, if sphere  A= pi*Df^2/4

    extern const double x_new_1; // For bigger domain build
    extern const double x_new_2;
    extern const double y_new_1;
    extern const double y_new_2;

    // Multiresolution parameters
    extern const double D_0;   // maximum core size
    extern const double dcmin; // minimum equivalent distance

    // Spatial Hashing Parameters
    extern const double cell_size; //size of cell
    extern const double xmin; 
    extern const double ymin;
    extern const double zmin;
    extern const double xmax;
    extern const double ymax;
    extern const double zmax;
    
    // DC PSE operators parameter ===================================== //
    extern const double beta;
    extern const double epsilon;
    extern const int l_2_0_2; // moment condition for Q(2,0) or Q(0,2) with 2nd order accuracy
    extern const int l_2_0_4; // moment condition for Q(2,0) or Q(0,2) with 4th order accuracy
    extern const int l_1_0_2; // moment condition for Q(1,0) or Q(0,1) with 2nd order accuracy
    extern const int l_1_0_4; // moment condition for Q(1,0) or Q(0,1) with 4th order accuracy
    extern const int l_0_0_3; // moment condition for Q(0,0) with 3rd order accuracy

    // ====== BOX Generation for post-processing extra data =========== //
    extern const int opt_gdata_vel; // =1 Using Biot-savart calculate Velocity on grid via Vorticity Strength on grid, =2 Remeshing Velocity from particle velocity

    // parameter for link list ======================================== //
    extern const int dim; // dimension
    extern const int skf; // Gauss spline
    extern int opt_neighbor;
    extern const int r_scale;

    // Tree option
    extern const double expTree;  // The percentage of tree root cell size expansion

    // FMM option
    extern const int P_max;          // Number of the expansion
    extern const int level_max;       // The tree level of the finest cell
    extern const int level_min;       // The tree level of the coarser cell
    extern const int n_max;          // Maximum number of the particle inside the cell
    extern const double tol_dist;   // The maximum tolerance of cell size

    // parameter for FMM ============================================== //
    extern const int lmax; // max of FMM levels
    extern const int npcm;
    extern const int nfwrd;
    extern const int nbmax;
    extern const int nbl1;
    extern const int nbl2;
    extern const int nbl3;
    extern const int nbl4;
    extern const int nbmrl;
    extern const double tol2;

    // FMM and FMMF option and their box size for solution ============ //
    extern const int ioptfmm;    // = 0 Direct, = 1 FMM
    extern const int n_s;        // 10 if using fmmo
    extern const int ndp;        // number of expansion
    extern const int icutoff;    // 0.singular ;  1.super (high-oder) algebraic  ; 2. Gaussian  ; 3 super Gaussian
    extern const int iopt_inter; //
    extern const int par_ext;    // =1(external flow case), =0(internal flow case)

    // REMESHING, DIFFUSION, SEARCHING, CONVECTION, STRETCHING, LES === //
    extern const int nrmsh;        // Remeshing every "nrmsh" timestep
    extern const int opt_remesh_W; // =1 using M4', =2 Using M6, =3 M4 isotropic (Chatelain,P [2005]), =4 P14
    extern const int par_split;    // 4,5,7,9 only, otherwise just spreading core size
    extern const int parsplit;     //
    extern const int xlimit;       // in order to eliminate far field splitting
    extern const int opt_sptadapt; // = 1 using spatial_adaptation for CSM, =0 withou spatial_adaptation
    extern const int opt_search;   // =1 using direct searching O(N^2); =2 using link_list O(N)
    extern const int it_start_les; // if  it_start_les =1, with zero gpx gpy gpz then velocity gradient dudx dudy ... = 0 at it =1, then strain(it =1)=0...=> Cr2 = NaN
    extern const int diffusion_opt;

    //VIBRATION PROPERTIES
    extern const int vib ;// =0 no vibration simulation, =1 vibration simulation
    extern double SpringConst;  // Spring constant
    extern double DamperConst; // Damper constant
    extern double mass; //Object mass
    extern double m_d; // Shen 2009, rho * (volume that occupied by body); for 2d i think its an area (?)
    extern double inertia;
    extern double tetha_dot;
    extern const double m_star;
    extern const double k_star;
    extern const double c_star; 
    extern const double t_star; 
    extern const double i_star;
    extern const double chi_star;
    extern const double U_star;
    extern const double tetha_nol;
    extern double gaya; //Gaya awal
    extern double tetha;
    extern double momen;


    // namespace geom
    // {
    // extern const int edge;
    // }
} // namespace Pars

#endif
