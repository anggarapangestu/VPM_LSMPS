#include "penalization.hpp"

// Velocity Penalization Main Manager
void penalization::get_penalization(Particle &p, const Body &b, int iter_step)
{
    /* Procedure:
       1. Collect all particle near the body (offset by Pars::numpen)
       2. Set up the penalization parameter (for initial iteration)
       3. Calculate the chi parameter
       4. Perform the vorticity penalization
       5. Update the penalization result into the original particle
    */
    
    // Internal variables
    Particle _p;             // Temporary particle for evaluation (Near body particle only)
    _p.num = 0;              // Initial particle number value
    std::vector<int> _idx;   // The original particle index list
    double *_xmin_pen = new double[2];
    double *_xmax_pen = new double[2];

    // PROCEDURE 1: Collect near body particle
    // ********************************************************************
    // Evaluate the extremes position for penalization evaluation
    _xmin_pen[0] = b.min_pos[0] - Pars::numpen * Pars::sigma;
    _xmin_pen[1] = b.min_pos[1] - Pars::numpen * Pars::sigma;
    _xmax_pen[0] = b.max_pos[0] + Pars::numpen * Pars::sigma;
    _xmax_pen[1] = b.max_pos[1] + Pars::numpen * Pars::sigma;

    /* ======== Determine Preliminary "minimum and maximum" grid position ========
       because maybe the body size is not fit to "body-penalization nodes/ except 6 
       additional nodes around pen-body" (biger or smaller)
       to make sure the 6 additional nodes will be outside of body we need pre-calculation:
        . . . . . . |     '     | . . . . .
        . . . . . . |     '     | . . . . .
        . . . . . . |     '     | . . . . .
        . . . . . . |     '     | . . . . .
                        center
    */

    // Penalization time manager
    clock_t t_pen = clock();
    printf("\nCalculating penalization ... \n");

    // Take only the particle inside the penalization evaluation domain
    for (size_t i = 0; i < p.num; i++)
    {
        if ((p.x[i] >= _xmin_pen[0]) && (p.x[i] <= _xmax_pen[0]) &&
            (p.y[i] >= _xmin_pen[1]) && (p.y[i] <= _xmax_pen[1]))
        {
            _idx.push_back(i);
            _p.num++;
            _p.x.push_back(p.x[i]);
            _p.y.push_back(p.y[i]);
            _p.u.push_back(p.u[i]);
            _p.v.push_back(p.v[i]);
            _p.s.push_back(p.s[i]);
            _p.gz.push_back(p.gz[i]);
            _p.vorticity.push_back(p.vorticity[i]);
        }
    }

    // PROCEDURE 2: Set up initial penalization parameter
    // ********************************************************************
    // Update lambda only at zero iteration step
    if (iter_step == 0){this->lambda_def();}

    // PROCEDURE 3: Perform the vorticity penalization
    // ********************************************************************
    // Calculate the chi
    this->kai_def(b, _p);

    // PROCEDURE 4: Set up initial penalization parameter
    // ********************************************************************
    // Perform the penalization
    if (Pars::iterative == 1){
        this->no_slip(_p, b);
    }
    // else if (Pars::iterative == 2){
    //     this->no_slip_iterative(_p, p, b, iter_step);
    // }

    // PROCEDURE 5: Update the penalization result into the original particle
    // ********************************************************************
    // Update the penalization value to original particle variable
    for (size_t i = 0; i < _p.num; i++)
    {
        p.u[_idx[i]] = _p.u[i];
        p.v[_idx[i]] = _p.v[i];
        p.gz[_idx[i]] = _p.gz[i];
        p.vorticity[_idx[i]] = _p.vorticity[i];
        p.chi[_idx[i]] = _p.chi[i];
    }

    // Display penalization computational time
    t_pen = clock() - t_pen;
    printf("<-> Penalization comp. time:           [%f s]\n", (double)t_pen/CLOCKS_PER_SEC);

    // Release heap variable
    delete[] _xmin_pen;
    delete[] _xmax_pen;
}
