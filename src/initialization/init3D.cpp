#include "initialization.hpp"

//  The real initialization function in 3D domain
void initialization::init_domain_3d(Particle &par)
{   
    // internal temporary variable
    double _x, _y, _z;
    int i, j, k, count;

    // number of particle in x and y direction
    int _nx = std::ceil(Pars::lxdom / Pars::sigma);     // Number of particle in x direction
    int _ny = std::ceil(Pars::lydom / Pars::sigma);     // Number of particle in y direction
    int _nz = std::ceil(Pars::lzdom / Pars::sigma);     // Number of particle in y direction
    int _nparticle = _nx * _ny * _nz;                   // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);
    par.z.resize(_nparticle, 0.0e0);

    // Generate the particle position
    count = 0;
    for (i = 0; i < _nx; i++)
    {
        for (j = 0; j < _ny; j++)
        {
            for (k = 0; k < _nz; k++){
                _x = -Pars::xdom + (i + 0.5)*Pars::sigma;
                _y = -Pars::sigma * _ny * 0.5 + (j + 0.5)*Pars::sigma;
                _z = -Pars::sigma * _nz * 0.5 + (k + 0.5)*Pars::sigma;
                par.x[count] = _x;
                par.y[count] = _y;
                par.z[count] = _z;
                count ++;
            }
        }
    }
    
    // Assign other particle properties
    par.s.resize(_nparticle, Pars::sigma);
    par.gz.resize(_nparticle, 0.0e0);
    par.u.resize(_nparticle, Pars::u_inf);
    par.v.resize(_nparticle, Pars::v_inf);
    par.w.resize(_nparticle, Pars::w_inf);
    par.isActive.resize(_nparticle, false);
    par.num = _nparticle;
}

/*
Note:
*/