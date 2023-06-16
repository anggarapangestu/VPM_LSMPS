#ifndef INCLUDED_GLOBAL
#define INCLUDED_GLOBAL

// #ifndef INCLUDED_UTILS
// #include "Utils.hpp"
// #endif

// Mains and basics
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>   // clock_t
#include <chrono>   // clock_t
#include <cmath>

// Additionals
#include <omp.h>    // pragma omp ...
#include <algorithm>
// #include <math.h>
// #include <complex>  // std::complex<>
// #include <numeric>

namespace Pars
{
// Predescription for each parameter variable
// [A] -> Already given in this program
// [I] -> Input for PARAMETER form user
// [O] -> Input for OPTION form user, if not will use default
// [C] -> Calculated form user input

    // ===================================================== //
    // +----------------- Simulation Flag -----------------+ //
    // ===================================================== //
    // Saving Flag
    extern const bool flag_save_body;       // The flag to save body data
    extern const bool flag_save_cell;       // The flag to save cell list data
    extern const bool flag_save_parameter;  // The flag to save simulation parameter
    extern const bool flag_save_sim_time;   // The flag to save simulation time and computational time
    extern const bool flag_save_residual;   // The flag for calculating and saving residual

    // Displaying Flag
    extern const bool flag_disp_stability;  // The flag to display stability parameter at each iteration
    extern const bool flag_disp_pred_time;  // The flag to display the prediction of total computational time until the end of simulation
    
    // Operator Flag
    extern const bool flag_pressure_calc;   // The flag to calculate pressure on post processing
    extern const bool flag_adaptive_dist;   // The flag for adaptive particle distribution

    // ================================================================= //
    // +----------------- Simulation Parameter Option -----------------+ //
    // ================================================================= //
    extern const int opt_sim_cont;  // [O] Option for starting simulation\
                                            0:= Start at initial;\
                                            1:= continue from the iteration ("cont_num") see the variable at "Computation Parameter"
    extern const int opt_init;      // [O] The initialization type option; \
                                            0:= Testing Resolution;\
                                            1:= Single Resolution;\
                                            2:= Multi-Res Single Block;\
                                            3:= Multi-Res Multi Block; // [Not Ready YET] \
                                            4:= Multi-Res Body Adjusted;
    extern const int opt_body;      // [O] The body type option; \
                                            1:= 2D - Circular cylinder, 3D - Sphere;\
                                            2:= Flat plate;\
                                            3:= Normal plate;\
                                            4:= Square cylinder;\
                                            5:= NACA airfoil;
    extern const int opt_neighbor;  // [O] The neighbor evaluation type option; \
                                            0:= Direct Neighbor Search;\
                                            1:= Linked List;\
                                            2:= Cell List;\
                                            3:= Spatial Hash;
    extern const int opt_pen_iter;  // [O] The brinkman iterative type option; \
                                            1 := classical brinkmann, \
                                            2 := iterative brinkmann

    // ============================================================== //
    // +----------------- Simulation Domain Option -----------------+ //
    // ============================================================== //
    extern const int DIM;        // [I] Spatial domain dimension; 2:= 2D simulation; 3:= 3D simulation
    extern const double lxdom;   // [I] Initial condition domain length of particle distribution on x-axis
    extern const double lydom;   // [I] Initial condition domain length of particle distribution on y-axis
    extern const double lzdom;   // [I] Initial domain length (z-axis)
    extern const double xdom;    // [I] Negative x-direction domain length 
    extern const double xcenter; // [I] Initial domain center (x-axis), default 0.0
    extern const double ycenter; // [I] Initial domain center (y-axis), default 0.0

    /*  Domain Illustration, see below.
                     |----------------------------------|
                     |                                  |
                     |                                  |
                     |     (xcenter, ycenter)           |
                     |<------------>*                   | lydom
                     |     xdom                         |
                     |                                  |
                     |                                  |
                     |----------------------------------|
                                      lxdom
    */
    
    // ===================================================================== //
    // +----------------- Incompressible Fluid Properties -----------------+ //
    // ===================================================================== //
    extern const double RE;     // [I] The reynold number (U)
    extern const double u_inf;  // [I] Freestream x-direction velocity
    extern const double v_inf;  // [I] Freestream y-direction velocity
    extern const double w_inf;  // [I] Freestream z-direction velocity
    extern const double RHO;    // [I] Fluid density
    extern const double U_inf;  // [C] Freestream velocity magnitude
    extern const double NU;     // [C] Kinematic viscousity
    extern const double MU;     // [C] Dynamic viscousity

    // =========================================================== //
    // ------------------ Computation Parameter ------------------ //
    // =========================================================== //
    // Basic parameter
    extern const double sigma;      // [I] Particle core size; ? default: 0.0025, dt~(1/2)*sigma - (1/3)*sigma
    extern const double dt;         // [I] Simulation time step; ? default:0.001, dt <= phi_s*sigma^2/vis (Ploumhans [2000]) OR dt = phi_s*sigma^2*Courant^2/vis, where 0 < Courant <= 1
    extern const double sim_time;   // [I] The total simulation time
    extern const int cont_num;      // [I] The continue data number ID

    // Stability Criteria
    extern const double phi_s;      // [A] (Ploumhans [2000]) for the Euler explicit scheme,phi_s = 0.595, !for the Adams–Bashforth 2 scheme, phi_s = 0.297, BUt can not use AB2 due to redistribution changes np--> need Biot-Savart again
    extern const double Courant;    // [C] Courant number(C): C = U_inf * dt / sigma, where 0 < C <= 1
    extern const double Diffusion;  // [C] Diffusion number(Phi): Phi =  NU * dt / sigma^2, where 0 < Phi <= 0.5
    // extern const double Courant =           // sigma = sqrt(dt*vis/phi_s)/Courant, where 0 < Courant <= 1
    //     std::sqrt(dt * vis / phi_s) / sigma;
    
    // Other utilities
    extern const int maxDigitLen;   // [C] Number of save name ID maximum digit
    extern const int nrmsh;         // [I] Step interval for remeshing evaluation

    // =========================================================== //
    // ------------------ Neighboring Parameter ------------------ //
    // =========================================================== //
    extern const double r_sup;      // [I] Neighbor evaluation support domain ratio size
    extern const double r_buff;     // [I] The buffer region radius factor
    extern const double body_ext;   // [I] The body extention distance to evaluate the chi and active particle
    extern const int max_level;     // [I] Number of resolution level
                                        // -> Level is counted from 0 as the larger particle size
                                        // -> Maximum level is the finest particle size
                                        // -> The number of resolution is 'max_level' + 1

    // ========================================================== //
    // ------------------ Solid Body Parameter ------------------ //
    // ========================================================== //
    // Whole Geometry Parameter
    extern const int n_a;           // [I] Number of body surface node
    extern const double ly;         // [I] flat plate: length, sphere: diameter (flat plate: length, sphere: diameter, Square: side.)
    extern const double lx;         // [I] sphere: diameter; chord*length ! lx*ly ! Flatplate, if sphere  A= pi*Df^2/4
    extern const double Df;         // [I] Reference length, flat plate: chord, sphere: diameter
    extern const double H_star;     // [I] Parameter for flat plate and normal plate: H*=H/L
    
    // Body velocity
    extern const int opt_u;         // [O] Body velocity option; 1= Suddenly start and stop (Instanteneosly appear), 2= smoothly start and stop, 3= moving
    extern const double ubody;      // [I] Body x-direction velocity
    extern const double vbody;      // [I] Body y-direction velocity

    // NACA Series parameter
    // NACA4-digit series: NACA(m*100)(p*10)(t*100) 2412,4612, Cylinder, or Flat plate
    extern const double m_a;        // [I] max thickness
    extern const double p_a;        // [I] maximum thickness position in percent of chord
    extern const double t_a;        // [I] is maximum thickness of airfoil in tenth of chord
    extern const double acx;        // [I] aerodynamic center location
    extern const double acy;        // [I] aerodynamic center location
    extern const double alpha_a;    // [I] alpha in degree , rotated angle (attached angle), ! Rotated airfoil

    // ============================================================ //
    // ------------------ Penalization Parameter ------------------ //
    // ============================================================ //
    extern const int opt_force_type;    // [O] The force calculation type;\
                                                    0:= No force calculation; \
                                                    1:= By penalization method; \
                                                    2:= By NOCA method
    extern const int opt_pen;           // [O] The penalization type option; \
                                                    1:= implicit penalization; (lambda = 1.0e4) \
                                                    2:= semi-implicit penalization,Rasmussen 2011 mixed with explicit; \
                                                    3:= Explicit scheme, Rasmussen 2011 (lambda = 1.0d0/dt)
    extern int opt_kaipen;              // [O] The χ (chi/kai) option; 0:= no recalculate kai; 1:= recalculated by kai~ ! to avoid singularities Rasmussen phD 2011
    extern double lambda;               // [O] Value of penalization
    extern const int numpen;            // [I] Number of points from surface for penalization domain;beberapa referensi ngasihnya 6
    extern const double hmollif;        // [I] The mollification length [2sqrt(2)*sigma] <- 2D, the half of it [sqrt(2)*sigma]

    // =================================================== //
    // ------------------ FMM Parameter ------------------ //
    // =================================================== //
    // Tree option
    extern const double expTree;    // [I] The percentage of tree root cell size expansion

    // FMM option
    extern const int ioptfmm;       // [O] Selection of the velocity calculation method;\
                                        0:= Direct calculation,\
                                        1:= FMM accelerated
    extern const int P_max;         // [O] Number of the expansion (P_max=10 for error in order of 10^-3)
    extern const int level_max;     // [I] The tree level of the finest cell
    extern const int level_min;     // [I] The tree level of the coarse cell
    extern const int n_max;         // [I] Maximum number of the particle inside the cell
    extern const double tol_dist;   // [A] The maximum tolerance of cell size (In evaluating adjacent cell)

    // Old FMM parameter
    extern const int lmax;          // [I] Max of FMM levels
    extern const int npcm;          // [I] Number of pole expansion
    extern const int nfwrd;         // std::pow(4,(lmax+1));   //nfwrd * (std::pow(4, (lmax + 1)) - 4) / 3;
    extern const int nbmax;         // Geometric Series formula
    extern const int nbmrl;         //nfwrd * 4096; //!!need to change if so that nbmrl >= 4^lmax
    extern const int nbl1;          // [A] Number of list 1 neighbor
    extern const int nbl2;          // [A] Number of list 2 neighbor
    extern const int nbl3;          // [A] Number of list 3 neighbor
    extern const int nbl4;          // [A] Number of list 4 neighbor
    extern const double tol2;       // [A] The maximum tolerance of cell size (In evaluating adjacent cell)
    extern const int n_s;           // [I] Number of maximum total particle inside a cell
    extern const int ndp;           // [I] The order of multipole expansion
    extern const int icutoff;       // [O] The redistribution function for direct potential calculation\
                                        0:= Singular,\
                                        1:= Super (high-oder) algebraic,\
                                        2:= Gaussian,\
                                        3:= Super Gaussian
                                    // Note:\
                                        -> if PSE calculation : take 2 or 3;\
                                        -> if Direct caluculation : take 1;\
                                        -> if FMM calculation : take 2 or 3

    // =========================================================== //
    // ------------------ Data Saving Parameter ------------------ //
    // =========================================================== //
    extern const double comtime_sf;     // % Save file frequently after how many [second]=[hour]*60*60, only in case of running is stopped suddenly, but nt_sf take long days
    extern const int opt_extra_data;    // 1:= not print extra data, print only common data
                                            // 2:= print common data and extra data at the same time,
                                            // 3:= print extra data from common data, MUST run case 1/2 first
    extern const int nt;                // [C] Total iterations number
    extern const int nt_data;           // [I] Total of data saved
    extern const int step_inv;          //std::ceil(nt / nt_data);   // [C] Step interval for saving data
    extern const int nt_sf;             // save data per ( nt_sf ) time step, for saving storage MEMORY and for case of running is stopped suddenly

    // ==================================================================== //
    // ------------------ Structural Vibration Parameter ------------------ //
    // ==================================================================== //
    // VIBRATION / STRUCTURAL PROPERTIES : !!STILL IN PROGRESS, RESULTS NOT GOOD ENOUGH!!
    extern const int vib;           // =0 no vibration simulation, =1 circular cylinder simulation, =2 flow induced pitch oscillations
    extern double SpringConst;      // Spring constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
    extern double DamperConst;      // Damper constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
    extern double mass;             //Object mass :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
    extern double inertia;          // Inertia :: !!STILL UNUSED!!
    //extern double m_d;            // !! rho times volume enclosed !! 
    extern double m_d;              // !! rho times volume enclosed !! UNUSED, SET TO 0.0 
    extern double tetha;            //!! ROTATION PARAMETER, STILL BUGGED
    extern double tetha_dot;        //!! ROTATION PARAMETER, STILL BUGGED
    extern const double m_star;     // NONDIMENSIONAL MASS
    extern const double k_star;     // NONDIMENSIONAL SPRING CONST
    extern const double c_star;     //NONDIMENSIONAL DAMPER CONST
    extern const double t_star;     //NONDIMENSIONAL TIME VARIABLE
    extern double gaya;             // !! FORCE INITIALIZATION !!
    extern double momen;            // !! MOMENT INITIALIZATION !!
    extern const double chi_star;   //!! ROTATION PARAMETER, STILL BUGGED
    extern const double i_star;     //!! ROTATION PARAMETER, STILL BUGGED
    extern const double U_star;     //!! ROTATION PARAMETER, STILL BUGGED
    extern const double tetha_nol;  //equilibrium angular position of spring [deg] //!! ROTATION PARAMETER, STILL BUGGED

} // namespace Pars

#endif
