#include "penalization.hpp"

// ================= No-slip BCs using PENALIZATION ======================
// TODO: Including 05 steps:
// *1st step: create_gpen.f90 Create Penalization grid nodes,
// *              Could be outside the main loop if non-deformable body.
// *2nd step: kai_def_gpen.f90 Define the Kai(chi) on penalization grid nodes,
// *              Could be outside the main loop with 1st step if non-deformable body.
// *3rd step: trans_vel_par2gpen.f90  Transfer particles velocity to penalization-grid-nodes
// *              velocity and vorticity (if using biot-savart: no need to transfer velocity )
// *4nd step: no_slip_bc_pen_gpen.f90 do penalization technique on penalization-grid
// *             from velocity on grid-pen, we calculate increment of velocity and vorticity penalization
// *5nd step: trans_incre_velvor_gpen2par.f90  tranfer increment of velocity and vorticity penalization on grid-pen to particles
// ======================================================

// ===================== STEP 02 ==============================
// TODO: 2nd step: kai_def_gpen.f90 penalization-grid defination: the Chi (kai), Could be outside the main loop with 1st step if fixed body.
// * - to calculate the indexes of each penalization-grid-node, to know it s inside or outside body and how far distance
// * outputs is _gpen variables and :h_min,n_gpen,ngx,ngy
// ============================================================

void penalization::kai_def(const Body &b, Particle &_p)
{
    // Calculate the nearest particle distance to body surface
    this->d_geom.distance_calc(_p,b);
    
    // Adjust the mollification parameter
    printf("<+> Calculate chi parameter\n");
    _p.chi.resize(_p.num, 0.0e0);
    for (int i = 0; i < _p.num; i++){
        // Chi calculation and evaluation
        if (_p.R[i] > Pars::hmollif) // Outside the body and mollification region
        {
            _p.chi[i] = 0.0e0;
        }
        else if (std::abs(_p.R[i]) <= Pars::hmollif) // Inside the mollification region (transition region)
        {
            _p.chi[i] = -0.5e0 * (-1.0e0 + _p.R[i] / Pars::hmollif + (1.0e0 / M_PI) * sin(M_PI * _p.R[i] / Pars::hmollif));
        }
        else if (_p.R[i] < -Pars::hmollif) // Inside the body region only
        {
            _p.chi[i] = 1.0e0;
        }
    }
    
    // Special mask function for semi-implicit scheme. [Rasmussen 2011]
    if (Pars::opt_kaipen == 1)
    {
        // Recalculating chi
        for (int i = 0; i < _p.x.size(); i++)
        {
            double kai_tilde;
            
            // Recalculation kai for semi-implicit brinkmann penalization (Combine the implicit and explicit)
            if (_p.chi[i] < 1.0e0){
                kai_tilde = _p.chi[i] / (Pars::dt * lambda * (1.0e0 - _p.chi[i]));
            }else{
                kai_tilde = _p.chi[i];
            }

            // Minimum between 1.0 and kai_tilde
            _p.chi[i] = kai_tilde < 1.0 ? kai_tilde : 1.0;
        }
    }
}

void penalization::lambda_def(){
    // TODO: Calculate penalization velocity u~, u^, u_lambda; 
    
    // Update the value of lambda; References [Rasmussen 2011]
    // IMPLICIT Scheme
    if (Pars::opt_pen == 1)
    {
        lambda = Pars::lambda;
    }
    // SEMI-IMPLICIT Scheme
    else if (Pars::opt_pen == 2)
    {
        Pars::opt_kaipen = 1;
        lambda = Pars::lambda;
    }
    // EXPLICIT Scheme
    else if (Pars::opt_pen == 3){
        lambda = 1.0e0 / Pars::dt;
    }

    // Update the value of Pars lambda
    Pars::lambda = lambda;
}

