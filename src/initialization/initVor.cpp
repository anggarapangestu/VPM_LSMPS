#include "initialization.hpp"


// ===================================================== //
// +--------------- Vorticity Initials ----------------+ //
// ===================================================== //
// #pragma region VORTICITY_INITIALS

/**
 *  @brief  Generate perlman vorticity as initial value.
 *         
 *  @param  _particle  Particle data for vorticity calculation.
 * 
 *  @return No return.
*/
void initialization::perlman_vorticity(Particle &par){
    /** Perlman Vorticity Field Formula (see. Rasmussen [2011])
     *            _
     *           /  (1-r^2)^7  ; r <= 1
     *   w(r) = <
     *           \_   0        ; r > 1
     * 
     *  The velocity solution to Perlman Vorticity
     *  >> x velocity  _
     *               /   -y/(16*r^2) * (1 - (1-r^2)^8)  ; r <= 1
     *     u(x,y) = <
     *               \_  -y/(16*r^2)                    ; r > 1
     * 
     *  >> y velocity  _
     *               /   x/(16*r^2) * (1 - (1-r^2)^8)  ; r <= 1
     *     v(x,y) = <
     *               \_  x/(16*r^2)                    ; r > 1
    */

    // Reserve the vorticity container
    par.vorticity.resize(par.num, 0.0e0);
    par.gz.resize(par.num, 0.0e0);

    // Calculate the vorticity at each particle
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Calculate the polar coordinate of radius
        double r2 = 0;
        r2 += (par.x[i]*par.x[i]);  // Component x
        r2 += (par.y[i]*par.y[i]);  // Component y
        #if (DIM == 3)
        r2 += (par.z[i]*par.z[i]);  // Component z
        #endif

        // Calculate the perlman vorticity
        if (r2 <= 1){
            par.vorticity[i] = std::pow((1.0 - r2), 7.0);
            par.gz[i] = par.vorticity[i] * par.s[i]*par.s[i];
        }else{
            par.vorticity[i] = 0;
            par.gz[i] = 0.0;
        }
    }
    return;
}


/**
 *  @brief  Calculate the analytical solution of perlman vorticity.
 *         
 *  @param  _particle  Particle data for velocity calculation.
 *  @param  _type   Choice to calculate the analytical solution:
 *                    > 0:= Calculate all analytical properties, 
 *                    > 1:= Velocity, 
 *                    > 2:= Stream function,
 *                    > 3:= Velocity differential
 * 
 *  @return No return.
*/
void initialization::perlman_velocity_solution(Particle &par, int type){
    /** Perlman Vorticity Field Formula
     *            _
     *           /  (1-r^2)^7  ; r <= 1
     *   w(r) = <
     *           \_   0        ; r > 1
     * 
     *  The velocity solution to Perlman Vorticity
     *  >> x velocity  _
     *               /   -y/(16*r^2) * (1 - (1-r^2)^8)  ; r <= 1
     *     u(x,y) = <
     *               \_  -y/(16*r^2)                    ; r > 1
     * 
     *  >> y velocity  _
     *               /   x/(16*r^2) * (1 - (1-r^2)^8)  ; r <= 1
     *     v(x,y) = <
     *               \_  x/(16*r^2)                    ; r > 1
     * 
     *  The perlman stream function field formula
     *              _
     *             /  (-4r^2 + 7r^4 - (28/3)r^6 + (35/4)r^8 - (28/5)r^10  ; r <= 1
     *   phi(r) = <      + (7/3)r^12 - (4/7)r^14 + (1/16)r^16) / 16
     *             \_   -(ln(r)+1.3589285714)/16    ; r > 1
     * 
     * 
     *  The velocity differential to Perlman Vorticity
     *  >> x velocity    _
     *                  /   xy/(8*r^4) * (1 - (1+7r^2)*(1-r^2)^7)  ; r <= 1
     *     dudx(x,y) = <
     *                  \_  xy/(8*r^4)                             ; r > 1
     *                   _
     *                  /   -1/(16*r^4) * (r^2-2y^2 + (14y^2r^2-r^2+r^4+2y^2)(1-r^2)^7)  ; r <= 1
     *     dudy(x,y) = <
     *                  \_  (y^2-x^2)/(16*r^4)                    ; r > 1
     * 
     *  >> y velocity  _
     *                  /   1/(16*r^4) * (r^2-2x^2 + (14x^2r^2-r^2+r^4+2x^2)(1-r^2)^7)  ; r <= 1
     *     dvdy(x,y) = <
     *                  \_  (y^2-x^2)/(16*r^4)                     ; r > 1
     *                   _
     *                  /   -xy/(8*r^4) * (1 - (1+7r^2)*(1-r^2)^7)  ; r <= 1
     *     dvdy(x,y) = <
     *                  \_  -xy/(8*r^4)                    ; r > 1
    */

    // Analytical solution data
    if (type == 0 || type == 1){
        par.u_a = std::vector<double> (par.num, 0.0e0);     // Analytical velocity in x direction
        par.v_a = std::vector<double> (par.num, 0.0e0);     // Analytical velocity in y direction
    }
    
    // Dirichlet boundary condition
    if (type == 0 || type == 2){
        par.phi_a = std::vector<double> (par.num, 0.0e0);   // Analytical stream function
    }

    // Neuman Boundary condition
    if (type == 0 || type == 3){
        par.dudx = std::vector<double> (par.num, 0.0e0);    // Differential in x for x velocity
        par.dudy = std::vector<double> (par.num, 0.0e0);    // Differential in y for x velocity
        par.dvdx = std::vector<double> (par.num, 0.0e0);    // Differential in x for y velocity
        par.dvdy = std::vector<double> (par.num, 0.0e0);    // Differential in y for y velocity
    }

    // Calculate the velocity solution at each particle
    if (type == 0){
        #pragma omp parallel for
        for (int i = 0; i < par.num; i++){
            // Calculate the polar coordinate of radius
            double r2 = 0;
            r2 += (par.x[i]*par.x[i]);  // Component x
            r2 += (par.y[i]*par.y[i]);  // Component y

            // Calculate the perlman velocity
            if (r2 <= 1){
                par.u_a[i] = -par.y[i] / (16.0*r2) * (1 - std::pow((1.0 - r2), 8.0));
                par.v_a[i] =  par.x[i] / (16.0*r2) * (1 - std::pow((1.0 - r2), 8.0));
            }else{
                par.u_a[i] = -par.y[i] / (16.0*r2);
                par.v_a[i] =  par.x[i] / (16.0*r2);
            }

            // Calculate the perlman stream function
            if (r2 <= 1){
                // double r4  = std::pow(r2,2);
                // double r6  = std::pow(r2,3);
                // double r8  = std::pow(r2,4);
                // double r10 = std::pow(r2,5);
                // double r12 = std::pow(r2,6);
                // double r14 = std::pow(r2,7);
                // double r16 = std::pow(r2,8);
                par.phi_a[i] = (1.0/16.0) * (- 4*r2 + 7*std::pow(r2,2)
                                    - ( 28.0 / 3.0  ) * std::pow(r2,3)
                                    + ( 35.0 / 4.0  ) * std::pow(r2,4)
                                    - ( 28.0 / 5.0  ) * std::pow(r2,5)
                                    + ( 7.0  / 3.0  ) * std::pow(r2,6)
                                    - ( 4.0  / 7.0  ) * std::pow(r2,7)
                                    + ( 1.0  / 16.0 ) * std::pow(r2,8)
                                    );
            }else{
                par.phi_a[i] = -(1.3589285714 + (std::log(r2)/2.0)) / 16.0;
            }

            // Calculate the perlman velocity differential
            if (r2 <= 1){
                // Alias the coordinate
                const double &x = par.x[i];
                const double &y = par.y[i];
                const double r4 = r2*r2;

                par.dudx[i] =  x*y/(8.0*r4) * (1 - (1 + 7*r2)*std::pow((1.0 - r2), 7.0));
                par.dudy[i] = -1.0/(16.0*r4) * (r2 - 2*y*y + (14*y*y*r2 - r2 + r4 + 2*y*y)*std::pow((1.0 - r2), 7.0));
                par.dvdx[i] =  1.0/(16.0*r4) * (r2 - 2*x*x + (14*x*x*r2 - r2 + r4 + 2*x*x)*std::pow((1.0 - r2), 7.0));
                par.dvdy[i] = -par.dudx[i];
            }else{
                // Alias the coordinate
                const double &x = par.x[i];
                const double &y = par.y[i];
                const double r4 = r2*r2;

                par.dudx[i] =  x*y/(8.0*r4);
                par.dudy[i] = (y*y - x*x)/(16.0*r4);
                par.dvdx[i] =  par.dudy[i];
                par.dvdy[i] = -par.dudx[i];
            }
        }
    }
    
    else if (type == 1){
        #pragma omp parallel for
        for (int i = 0; i < par.num; i++){
            // Calculate the polar coordinate of radius
            double r2 = 0;
            r2 += (par.x[i]*par.x[i]);  // Component x
            r2 += (par.y[i]*par.y[i]);  // Component y

            // Calculate the perlman velocity
            if (r2 <= 1){
                par.u_a[i] = -par.y[i] / (16.0*r2) * (1 - std::pow((1.0 - r2), 8.0));
                par.v_a[i] =  par.x[i] / (16.0*r2) * (1 - std::pow((1.0 - r2), 8.0));
            }else{
                par.u_a[i] = -par.y[i] / (16.0*r2);
                par.v_a[i] =  par.x[i] / (16.0*r2);
            }
        }
    }

    else if (type == 2){
        #pragma omp parallel for
        for (int i = 0; i < par.num; i++){
            // Calculate the polar coordinate of radius
            double r2 = 0;
            r2 += (par.x[i]*par.x[i]);  // Component x
            r2 += (par.y[i]*par.y[i]);  // Component y

            // Calculate the perlman stream function
            if (r2 <= 1){
                par.phi_a[i] = (1.0/16.0) * (- 4*r2 + 7*std::pow(r2,2)
                                    - ( 28.0 / 3.0  ) * std::pow(r2,3)
                                    + ( 35.0 / 4.0  ) * std::pow(r2,4)
                                    - ( 28.0 / 5.0  ) * std::pow(r2,5)
                                    + ( 7.0  / 3.0  ) * std::pow(r2,6)
                                    - ( 4.0  / 7.0  ) * std::pow(r2,7)
                                    + ( 1.0  / 16.0 ) * std::pow(r2,8)
                                    );
            }else{
                par.phi_a[i] = -(1.3589285714 + (std::log(r2)/2.0)) / 16.0;
            }
        }
    }

    if (type == 3){
        #pragma omp parallel for
        for (int i = 0; i < par.num; i++){
            // Calculate the polar coordinate of radius
            double r2 = 0;
            r2 += (par.x[i]*par.x[i]);  // Component x
            r2 += (par.y[i]*par.y[i]);  // Component y

            // Calculate the perlman velocity differential
            if (r2 <= 1){
                // Alias the coordinate
                const double &x = par.x[i];
                const double &y = par.y[i];
                const double r4 = r2*r2;

                par.dudx[i] =  x*y/(8.0*r4) * (1 - (1 + 7*r2)*std::pow((1.0 - r2), 7.0));
                par.dudy[i] = -1.0/(16.0*r4) * (r2 - 2*y*y + (14*y*y*r2 - r2 + r4 + 2*y*y)*std::pow((1.0 - r2), 7.0));
                par.dvdx[i] =  1.0/(16.0*r4) * (r2 - 2*x*x + (14*x*x*r2 - r2 + r4 + 2*x*x)*std::pow((1.0 - r2), 7.0));
                par.dvdy[i] = -par.dudx[i];
            }else{
                // Alias the coordinate
                const double &x = par.x[i];
                const double &y = par.y[i];
                const double r4 = r2*r2;

                par.dudx[i] =  x*y/(8.0*r4);
                par.dudy[i] = (y*y - x*x)/(16.0*r4);
                par.dvdx[i] =  par.dudy[i];
                par.dvdy[i] = -par.dudx[i];
            }
        }
    }
    return;
}


/**
 *  @brief  Calculate the eliptical vorticity as initial value.
 *         
 *  @param  _particle  Particle data for velocity calculation.
 * 
 *  @return No return.
*/
void initialization::eliptic_vorticity(Particle &par){
    /** The eliptic_vorticity (see. Rossineli [2015])
     *             _
     *            /   Vm(1 - f(r/R0))  ; r/R0 <= 1
     *    w(r) = <
     *            \_         0         ; r > 1
     *  where:
     *    f(z) = exp( -(q0/z) * exp(1/(z-1)) )
    */
    
    // Internal parameter
    double R0 = 0.8;        // Base radius
    double Vm = 20.0;       // Maximum vorticity
    double q0 = 2.56085;    // Smoothing constant

    // Reserve the vorticity data
    par.vorticity.resize(par.num);
    par.gz.resize(par.num);

    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Calculate the polar coordinate of radius
        double r2 = 0;
        r2 += 1*(par.x[i]*par.x[i]);    // Component x
        r2 += 0.25*(par.y[i]*par.y[i]);      // Component y

        // Calculate the other parameter
        double r = std::sqrt(r2);       // The elliptical radius
        double z = r/R0;               // The radius ratio

        if (z == 0.0){
            // Singular value at z = 0.0
            par.vorticity[i] = Vm;
            par.gz[i] = Vm * par.s[i]*par.s[i];
        }
        else if (z < 1.0){
            double f_q = std::exp(-(q0/z) * std::exp(1/(z-1)));
            par.vorticity[i] = Vm * (1 - f_q);
            par.gz[i] = par.vorticity[i] * par.s[i]*par.s[i];
        }else{
            par.vorticity[i] = 0.0;
            par.gz[i] = 0.0;
        }
    }
    return;
}


/**
 *  @brief  Generate lamb oseen vorticity as initial value. (based on the Rossi (2015))
 *         
 *  @param  _particle  Particle data for vorticity calculation.
 * 
 *  @return No return.
*/
void initialization::lamb_oseen_vorticity(Particle &par){
    /** Lamb Oseen Vorticity Field Formula (see. Rossi [2015])
     * 
     *  Basis value
     *   Set value : Re (Reynolds Number), 
     *               L (Vortex scaling length Half-Width of gaussian), 
     *               nu (Flow viscousity)
     *   
     *  Basis Parameter
     *  [1]      Re * nu
     *   W_0 =  ----------  -> Peak vorticity @t=0
     *           pi * L^2
     * 
     *  [2]            Re
     *   f(t) = -----------------   -> Time parameter function
     *           4*pi*t*W_0 + Re
     *  [3]
     *   g(r,t) = (r/L)^2 * f(t)    -> Spatial & time paremeter function
     * 
     * 
     *  Time variating properties (analytical solution)
     *   w(r,t) = W_0 * f(t) * exp(-g(r,t))
     *   u(r,t) = W_0 * (L^2 / (2*r)) * [1 - exp(-g(r,t))]
     * 
     *  Initial properties
     *   w(r,0) = W_0 * exp(-(r/L)^2)
     *   u(r,0) = W_0 * (L^2 / (2*r)) * [1 - exp(-(r/L)^2)]
     * 
     * Setting for Maximum domain size (IMPORTANT for simulation parameter setting !!!)
     *   a = Tolerance of dominant vorticity (here set to be 10^-6)
     *   T = Final simulation time
     * 
     *   R_min = sqrt(-(ln(a)) * L^2 / f(T))
     * And domain size must be x_dom > 2*R_min
    */

    // Basic Parameter
    const double L = 0.2;
    const double W_0 = Pars::RE * Pars::NU / (M_PI * L * L);

    // For initial configuration
    //  RE: 100 (higher cause higher peak), 
    //  Nu: 0.005 (higher cause higher peak)
    //  L: 0.2 (Need for good adjustment)
    // Then we get at T=20 the min R is 2.5 -> Domain size must be 5x5

    // Reserve the vorticity container
    par.vorticity.resize(par.num, 0.0e0);
    par.gz.resize(par.num, 0.0e0);

    // Calculate the vorticity at each particle
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Calculate the polar coordinate of radius
        double r2 = 0;
        r2 += (par.x[i]*par.x[i]);  // Component x
        r2 += (par.y[i]*par.y[i]);  // Component y

        // Calculate the lamb oseen vorticity [w(r,0) = W_0 * exp(-(r/L)^2)]
        par.vorticity[i] = W_0 * (std::exp(-r2 / (L * L)));
        par.gz[i] = par.vorticity[i] * par.s[i]*par.s[i];
    }
    return;
}

/**
 *  @brief  Calculation of the lamb oseen vorticity analytical.
 *         
 *  @param  _particle  Particle data for vorticity calculation.
 *  @param  _time  Current simulation time step.
 *  @param  _type  Type of the analytical data calculation.
 *                 [1:= Calculate all parameter]
 * 
 *  @return No return.
*/
void initialization::lamb_oseen_solution(Particle &par, double _t, int _type){
    // For initial configuration
    //  RE: 100 (higher cause higher peak), 
    //  Nu: 0.005 (higher cause higher peak)
    //  L: 0.2 (Need for good adjustment)
    // Then we get at T=20 the min R is 2.5 -> Domain size must be 5x5
    
    // Basic Parameter
    // const double L = 0.2;   // Make sure this value similar to the one in vorticity generation
    const double L = Pars::L_LO;
    const double W_0 = Pars::RE * Pars::NU / (M_PI * L * L);

    // Function Parameter
    const double f_t = Pars::RE / (4*M_PI*_t*W_0 + Pars::RE);

    // Reserve the vorticity container
    par.vort_a.clear(); par.vort_a.resize(par.num, 0.0e0);
    par.u_a.clear(); par.u_a.resize(par.num, 0.0e0);
    par.v_a.clear(); par.v_a.resize(par.num, 0.0e0);

    // Calculate the vorticity at each particle
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Calculate the polar coordinate of radius
        double r2 = 0;
        r2 += (par.x[i] * par.x[i]);  // Component x
        r2 += (par.y[i] * par.y[i]);  // Component y
        double _r = sqrt(r2);

        // Function parameter g
        double g_r_t = f_t * r2 / (L * L);

        // Orientation of particle motion
        double cosT = -par.y[i] / _r;
        double sinT = par.x[i] / _r;

        // Calculate the lamb oseen vorticity [w(r,0) = W_0 * exp(-(r/L)^2)]
        par.vort_a[i] = W_0 * f_t * std::exp(-g_r_t);
        double Vel = W_0 * L * L / (2*_r) * (1 - std::exp(-g_r_t));

        par.u_a[i] = Vel * cosT;
        par.v_a[i] = Vel * sinT;
    }
    
    return;
}


// The real ring vortex simulation
// /**
//  *  @brief  Calculate the vortex ring initial vorticity.
//  *         
//  *  @param  _particle  Particle data for velocity calculation.
//  * 
//  *  @return No return.
// */
// void initialization::ring_vortex(Particle &par){
//     /** The ring vortex initials
//      *   vorticity distribution
//      *    w(r) = G / (2*pi*sig^2) * exp(-r^2/(2*sig^2))
//      * 
//      *   where:
//      *    sig = 0.05 * R
//      * 
//      *  It gives a ring vortex with inner radius of 0.25 R
//      *  NOTE: Particle size need smaller than sigma
//     */

//     // // Number of vortex ring
//     // int NUM = 2;
    
//     // Reserve the vorticity data
//     par.vortx.resize(par.num);
//     par.vorty.resize(par.num);
//     par.vortz.resize(par.num);
//     par.vorticity.resize(par.num);

//     // **Geometrical Parameter
//     double Rout = 1.0;          // Outer radius (the ring circle radius)
//     double Rin = 0.25*Rout;     // Inner radius (the ring cross section radius)
//     std::vector<double> xCen = {0.0,0.0,0.0};       // The center location
//     std::vector<double> xNorm = {1.0,1.0,1.0};      // The normal vector
//     // std::vector<double> xCen1 = {0.0, Rout/2, Rout/8};       // The center location
//     // std::vector<double> xCen2 = {0.0, -Rout/2, -Rout/8};       // The center location
//     // std::vector<double> xNorm1 = {0.0,0.0,1.0};      // The normal vector
//     // std::vector<double> xNorm2 = {0.0,0.0,-1.0};      // The normal vector
//     {// Normalize the vector of ring normal direction
//         double NORM = 0.0;
//         for (auto &elm : xNorm){
//             NORM += elm*elm;
//         }
//         NORM = std::sqrt(NORM);
//         for (auto &elm : xNorm){
//             elm = elm/NORM;
//         }
//     }

//     // **Vorticity parameter
//     double G = 1.0;         // Ring vorticiy scaling
//     double R = Rout;        // Base ring radius
//     double sig = 0.05 * R;  // Smoothing constant

//     // **Ring Vortex evaluation
//     #pragma omp parallel for
//     for (int i = 0; i < par.num; i++){
//         // Calculate the distance vector to ring center
//         std::vector<double> _r = {0,0,0};
//         _r[0] = par.x[i] - xCen[0];
//         _r[1] = par.y[i] - xCen[1];
//         _r[2] = par.z[i] - xCen[2];

//         // Calculate the distance vector projection to ring plane and the normal distance
//         /** The illustration of distance vector triangle
//          *                          .
//          *             Dist   .     |  ^ Norm
//          *            O .___________|  |
//          *                  ---->
//          *                   Proj
//         */
//         double _dist = 0.0, _norm = 0.0, _proj = 0.0;
//         // The norm distance is the dot product to ring normal vector
//         for (int d = 0; d < 3; d++){
//             _dist += _r[d] * _r[d];
//             _norm += xNorm[d] * _r[d];
//         }
//         _dist = std::sqrt(_dist);
//         _proj = std::sqrt(_dist*_dist + _norm*_norm);


//         // Calculate the local radius
//         double r2 = 0;
//         r2 += std::pow(_proj-R, 2.0);  // Component x
//         r2 += std::pow(_norm, 2.0);    // Component y
//         double r = std::sqrt(r2);      // The local radius

//         // Calculate the tangential vector of current particle relative 
//         std::vector<double> _tgn = {0.0, 0.0, 0.0};
//         // The tangential vector is cross product between norm and distance vector
//         _tgn[0] = xNorm[1]*_r[2] - xNorm[2]*_r[1];
//         _tgn[1] = xNorm[2]*_r[0] - xNorm[0]*_r[2];
//         _tgn[2] = xNorm[0]*_r[1] - xNorm[1]*_r[0];
//         {// Normalize the tangential vector
//             double NORM = 0.0;
//             for (auto &elm : _tgn){
//                 NORM += elm*elm;
//             }
//             NORM = std::sqrt(NORM);
//             for (auto &elm : _tgn){
//                 elm = elm/NORM;
//             }
//         }

//         // Calculate the vorticity
//         if (r < Rin){
//             double omega = G / (2*M_PI*sig*sig) * std::exp(-(r*r)/(2*sig*sig));
//             par.vorticity[i] = omega;
//             par.vortx[i] = omega * _tgn[0];
//             par.vorty[i] = omega * _tgn[1];
//             par.vortz[i] = omega * _tgn[2];
//         }else{
//             par.vortx[i] = 0.0;
//             par.vorty[i] = 0.0;
//             par.vortz[i] = 0.0;
//             par.vorticity[i] = 0.0;
//         }
//     }
//     return;
// }

/**
 *  @brief  Calculate the vortex ring initial vorticity.
 *         
 *  @param  _particle  Particle data for velocity calculation.
 * 
 *  @return No return.
*/
void initialization::ring_vortex(Particle &par){
    /** The ring vortex initials
     *   vorticity distribution
     *    w(r) = G / (2*pi*sig^2) * exp(-r^2/(2*sig^2))
     * 
     *   where:
     *    sig = 0.05 * R
     * 
     *  It gives a ring vortex with inner radius of 0.25 R
     *  NOTE: Particle size need smaller than sigma
    */
    
    // Number of ring
    int NUM = 2;

    // Reserve the vorticity data
    par.vortx.resize(par.num);
    par.vorty.resize(par.num);
    par.vortz.resize(par.num);
    par.vorticity.resize(par.num);

    // **Geometrical Parameter
    std::vector<double> Ro(NUM);    // Outer radius (the ring circle radius)
    std::vector<double> Ri(NUM);    // Outer radius (the ring circle radius)
    std::vector<std::vector<double>> xCen(NUM);     // The center location
    std::vector<std::vector<double>> xNorm(NUM);    // The normal vector (determine the vector direction)
    
    // Parameter of First Ring
    Ro[0] = 1.0;                                // Outer radius can be modify
    Ri[0] = Ro[0] / 4.0;
    xCen[0]  = {0.0, Ro[0]/2.0, Ro[0]/8.0};     // Define the center position
    xNorm[0] = {0.0, 0.0, -1.0};                 // Define the normal direction

    // Parameter of Second Ring
    Ro[1] = 1.0;                                // Outer radius can be modify
    Ri[1] = Ro[1] / 4.0;
    xCen[1]  = {0.0, -Ro[1]/2.0, -Ro[1]/8.0};   // Define the center position
    xNorm[1] = {0.0, 0.0, 1.0};                // Define the normal direction
    
    // Normalize the vector of ring normal direction
    for (int i = 0; i < NUM; i++){
        double NORM = 0.0;
        for (auto &elm : xNorm[i]){
            NORM += elm*elm;
        }
        NORM = std::sqrt(NORM);
        for (auto &elm : xNorm[i]){
            elm = elm/NORM;
        }
    }

    // **Vorticity parameter
    double G = 1.0;         // Ring vorticiy scaling
    std::vector<double> R = Ro;  // Base ring radius
    std::vector<double> sig = R; // Smoothing constant
    for (int i = 0; i < NUM; i++){
        sig[i] *= 0.05;
    }

    // **Ring Vortex evaluation
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Main parameter
        double dist = Pars::lxdom;  // Initial distance
        int part = 0;               // Which ring part it correlated to
        std::vector<double> _R = {0.0, 0.0, 0.0};       // The distance vector

        // Check the distance toward each ring vortex
        for (int k = 0; k < NUM; k++){
            // Calculate the distance vector to ring center
            std::vector<double> _r = {0,0,0};
            _r[0] = par.x[i] - xCen[k][0];
            _r[1] = par.y[i] - xCen[k][1];
            _r[2] = par.z[i] - xCen[k][2];

            // Calculate the distance vector projection to ring plane and the normal distance
            /** The illustration of distance vector triangle
             *                          .
             *             Dist   .     |  ^ Norm
             *            O .___________|  |
             *                  ---->
             *                   Proj
            */
            double _dist = 0.0, _norm = 0.0, _proj = 0.0;
            // The norm distance is the dot product to ring normal vector
            for (int d = 0; d < 3; d++){
                _dist += _r[d] * _r[d];
                _norm += xNorm[k][d] * _r[d];
            }
            _dist = std::sqrt(_dist);
            _proj = std::sqrt(_dist*_dist + _norm*_norm);

            // Calculate the local radius
            double r2 = 0;
            r2 += std::pow(_proj - R[k], 2.0);  // Component x
            r2 += std::pow(_norm, 2.0);    // Component y
            double r = std::sqrt(r2);      // The local radius

            // Change the part
            if (r < dist){
                dist = r;
                part = k;
                _R = _r;
            }
        }
        
        // Calculate the tangential vector of current particle relative 
        std::vector<double> _tgn = {0.0, 0.0, 0.0};     // The tangential vector

        // The tangential vector is cross product between norm and distance vector
        _tgn[0] = xNorm[part][1]*_R[2] - xNorm[part][2]*_R[1];
        _tgn[1] = xNorm[part][2]*_R[0] - xNorm[part][0]*_R[2];
        _tgn[2] = xNorm[part][0]*_R[1] - xNorm[part][1]*_R[0];
        {// Normalize the tangential vector
            double NORM = 0.0;
            for (auto &elm : _tgn){
                NORM += elm*elm;
            }
            NORM = std::sqrt(NORM);
            for (auto &elm : _tgn){
                elm = elm/NORM;
            }
        }

        // Calculate the vorticity
        if (dist < Ri[part]){
            const double &s = sig[part];
            double omega = G / (2*M_PI*s*s) * std::exp(-(dist*dist)/(2*s*s));
            par.vorticity[i] = omega;
            par.vortx[i] = omega * _tgn[0];
            par.vorty[i] = omega * _tgn[1];
            par.vortz[i] = omega * _tgn[2];
        }else{
            par.vortx[i] = 0.0;
            par.vorty[i] = 0.0;
            par.vortz[i] = 0.0;
            par.vorticity[i] = 0.0;
        }
    }
    return;
}

/**
 *  @brief  Calculate the vortex ring initial vorticity.
 *         
 *  @param  _particle  Particle data for velocity calculation.
 * 
 *  @return No return.
*/
void initialization::ring_vortex_oblique(Particle &par){
    /** The ring vortex initials
     *   vorticity distribution
     *    w(r) = G / (2*pi*sig^2) * exp(-r^2/(2*sig^2))
     * 
     *   where:
     *    sig = 0.1 * R
     * 
     *  It gives a ring vortex with inner radius of 0.3 R
     *  NOTE: Particle size need smaller than sigma
    */
    
    // Number of ring
    int NUM = 2;

    // Reserve the vorticity data
    par.vortx.resize(par.num);
    par.vorty.resize(par.num);
    par.vortz.resize(par.num);
    par.vorticity.resize(par.num);

    // **Geometrical Parameter
    std::vector<double> Ro(NUM);    // Outer radius (the ring circle radius)
    std::vector<double> Ri(NUM);    // Outer radius (the ring circle radius)
    std::vector<std::vector<double>> xCen(NUM);     // The center location
    std::vector<std::vector<double>> xNorm(NUM);    // The normal vector (determine the vector direction)

    // Note
    // All vortex is facing down and tilted for 15 degree
    // The distance of center between 2 is 2.7R or 1.35 each

    double sinT = std::sin(15.0*M_PI/180.0);
    double cosT = std::cos(15.0*M_PI/180.0);
    
    // Parameter of First Ring (On left -> negative x)
    Ro[0] = 1.0;                                // Outer radius can be modify
    Ri[0] = Ro[0] / 2.5;  // The significant value of gaussian radius
    xCen[0]  = {-1.35, 0.0, 1.0};               // Define the center position
    xNorm[0] = {sinT, 0.0, -cosT};                 // Define the normal direction

    // Parameter of Second Ring (On right -> positive x)
    Ro[1] = 1.0;                                // Outer radius can be modify
    Ri[1] = Ro[1] / 2.5;
    xCen[1]  = {1.35, 0.0, 1.0};                // Define the center position
    xNorm[1] = {-sinT, 0.0, -cosT};                // Define the normal direction
    
    // Normalize the vector of ring normal direction
    for (int i = 0; i < NUM; i++){
        double NORM = 0.0;
        for (auto &elm : xNorm[i]){
            NORM += elm*elm;
        }
        NORM = std::sqrt(NORM);
        for (auto &elm : xNorm[i]){
            elm = elm/NORM;
        }
    }

    // **Vorticity parameter
    double G = 1.0;         // Ring vorticiy scaling
    std::vector<double> R = Ro;  // Base ring radius
    std::vector<double> sig = R; // Smoothing constant
    for (int i = 0; i < NUM; i++){
        sig[i] *= 0.1;
    }

    // **Ring Vortex evaluation
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Main parameter
        double dist = Pars::lxdom;  // Initial distance
        int part = 0;               // Which ring part it correlated to
        std::vector<double> _R = {0.0, 0.0, 0.0};       // The distance vector

        // Check the distance toward each ring vortex
        for (int k = 0; k < NUM; k++){
            // Calculate the distance vector to ring center
            std::vector<double> _r = {0,0,0};
            _r[0] = par.x[i] - xCen[k][0];
            _r[1] = par.y[i] - xCen[k][1];
            _r[2] = par.z[i] - xCen[k][2];

            // Calculate the distance vector projection to ring plane and the normal distance
            /** The illustration of distance vector triangle
             *                          .
             *             Dist   .     |  ^ Norm
             *            O .___________|  |
             *                  ---->
             *                   Proj
            */
            double _dist = 0.0, _norm = 0.0, _proj = 0.0;
            // The norm distance is the dot product to ring normal vector
            for (int d = 0; d < 3; d++){
                _dist += _r[d] * _r[d];
                _norm += xNorm[k][d] * _r[d];
            }
            _dist = std::sqrt(_dist);
            _proj = std::sqrt(_dist*_dist + _norm*_norm);

            // Calculate the local radius
            double r2 = 0;
            r2 += std::pow(_proj - R[k], 2.0);  // Component x
            r2 += std::pow(_norm, 2.0);    // Component y
            double r = std::sqrt(r2);      // The local radius

            // Change the part
            if (r < dist){
                dist = r;
                part = k;
                _R = _r;
            }
        }
        
        // Calculate the tangential vector of current particle relative 
        std::vector<double> _tgn = {0.0, 0.0, 0.0};     // The tangential vector

        // The tangential vector is cross product between norm and distance vector
        _tgn[0] = xNorm[part][1]*_R[2] - xNorm[part][2]*_R[1];
        _tgn[1] = xNorm[part][2]*_R[0] - xNorm[part][0]*_R[2];
        _tgn[2] = xNorm[part][0]*_R[1] - xNorm[part][1]*_R[0];
        {// Normalize the tangential vector
            double NORM = 0.0;
            for (auto &elm : _tgn){
                NORM += elm*elm;
            }
            NORM = std::sqrt(NORM);
            for (auto &elm : _tgn){
                elm = elm/NORM;
            }
        }

        // Calculate the vorticity
        if (dist < Ri[part]){
            const double &s = sig[part];
            double omega = G / (2*M_PI*s*s) * std::exp(-(dist*dist)/(2*s*s));
            par.vorticity[i] = omega;
            par.vortx[i] = omega * _tgn[0];
            par.vorty[i] = omega * _tgn[1];
            par.vortz[i] = omega * _tgn[2];
        }else{
            par.vortx[i] = 0.0;
            par.vorty[i] = 0.0;
            par.vortz[i] = 0.0;
            par.vorticity[i] = 0.0;
        }
    }
    return;
}

// #pragma endregion


// ==================================================
// +--------------- TESTING FUNCTION ---------------+
// ==================================================
/** All testing function store the data in 2D
 *  - Coordinate at x, y
 *  - Function in   phi_a (analytical)
 *  - Source in     F
*/

// GAUSSIAN FUNCTION
void initialization::test_0(Particle &par){
    // A gaussian function: phi(x,y) = exp (-(x^2 + y^2))   [Solution]
    // A gaussian function: p(x,y) = -4*(1-(x^2 + y^2)) * exp(-(x^2 + y^2))     [Source]
    // Reserve the container
    par.phi_a = std::vector<double>(par.num, 0.0e0);    // The analytical potential
    par.F = std::vector<double>(par.num, 0.0e0);        // The source of poisson

    par.gx = std::vector<double>(par.num, 0.0e0);       // Diferential toward x
    par.gy = std::vector<double>(par.num, 0.0e0);       // Diferential toward y

    // Calculate the function and the source
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Calculate the polar coordinate of radius
        double r2 = 0;
        r2 += (par.x[i]*par.x[i]);  // Component x
        r2 += (par.y[i]*par.y[i]);  // Component y

        // Calculate the perlman vorticity
        double phi = std::exp(-r2);
        par.phi_a[i] = phi;
        par.F[i] = -4 * (1 - r2) * phi;

        // For neuman boundary condition
        par.gx[i] = -2*par.x[i]*phi;
        par.gy[i] = -2*par.y[i]*phi;
    }
    return;
}

// POLYNOMIAL FUNCTION
void initialization::test_1(Particle &par){
    // A simple polynomial function: phi(x,y) = (x^2 * y^2) [Solution]
    // A simple polynomial function: f(x,y) = 2(x^2 + y^2)  [Source]
    // Reserve the container
    par.phi_a = std::vector<double>(par.num, 0.0e0);          // The potential
    par.F = std::vector<double>(par.num, 0.0e0);    // The source of poisson

    par.gx = std::vector<double>(par.num, 0.0e0);       // Diferential toward x
    par.gy = std::vector<double>(par.num, 0.0e0);       // Diferential toward y

    // Calculate the function and the source
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Alias the value
        const double &x = par.x[i];
        const double &y = par.y[i];

        // Calculate the perlman vorticity
        par.phi_a[i] = std::pow(x*y, 2);
        par.F[i] = 2 * (x*x + y*y);
        
        // For neuman boundary condition
        par.gx[i] = 2*x*y*y;
        par.gy[i] = 2*x*x*y;
    }
    return;
}

// MORE COMPLICATED FUNCTION
void initialization::test_2(Particle &par){
    // A trigonometric polynomial function: phi(x,y) = 2(x^2)sin^2(y) + 4(x^2)y + 3(y^2)cos(2x)     [Solution]
    // A trigonometric polynomial function: f(x,y) = 4sin^2(y) + 4(x^2)cos(2y) + 8y + (6-12y^2)cos(2x)  [Source]
    // Reserve the container
    par.phi_a = std::vector<double>(par.num, 0.0e0);          // The potential
    par.F = std::vector<double>(par.num, 0.0e0);    // The source of poisson

    par.gx = std::vector<double>(par.num, 0.0e0);       // Diferential toward x
    par.gy = std::vector<double>(par.num, 0.0e0);       // Diferential toward y

    // Calculate the function and the source
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Alias the value
        const double &x = par.x[i];
        const double &y = par.y[i];

        // Calculate the perlman vorticity
        par.phi_a[i] = 2*std::pow(x*sin(y), 2) + 4*x*x*y + 3*y*y*cos(2*x);
        par.F[i] = 4*std::pow(sin(y), 2) + 4*x*x*cos(2*y) + 8*y + (6-12*y*y)*cos(2*x);

        // For neuman boundary condition
        par.gx[i] = 4*x*std::pow(sin(y), 2) + 8*x*y - 6*y*y*sin(2*x);
        par.gy[i] = 2*x*x*sin(2*y) + 4*x*x + 6*y*cos(2*x);
    }
    return;
}

// Laplace solution
void initialization::laplace(Particle &par){
    // The basic solution through the laplace problem 
    // V(x,y) = 2*V_0/pi * arctan(sin(pi*(x/L)) / sinh (pi*(y/L)));

    // An internal variable for the laplace solution
    double L = Pars::lxdom;         // The domain side size
    double sx = 0.5*Pars::lxdom;        // The x boundary shift
    double sy = 0.5*Pars::lydom;        // The y boundary shift
    double V_0 = 10.0;       // The potential at the left side (x=-sx)
    int _P = 100;

    // Reserve the container
    par.phi_a = std::vector<double>(par.num, 0.0e0);
    par.F     = std::vector<double>(par.num, 0.0e0);
    
    // Calculate the analytical solution of potential
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Alias the value
        double x = par.x[i] + sx;
        double y = par.y[i] + sy;
        // double num = std::sin(pi*(y/L));
        // double den = std::sinh(pi*(x/L));
        // double V = 2*V_0/pi * std::atan(num / den);
        // if (std::isnan(V)){
        //     V = 2*V_0/pi * std::atan(1.0);
        // }

        double V = 0.0;
        for (int n = 1; n < _P; n++){
            V += (2*V_0/(n*M_PI))*((1-std::cos(M_PI *n))/std::sinh(n*M_PI))*
                 (std::sin(n*M_PI*y/L)*std::sinh(n*M_PI*x/L));
        } 

        // Calculate the analytical solution
        par.phi_a[i] = V;
    }
    return;
}