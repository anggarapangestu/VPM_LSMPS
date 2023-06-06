#include "neighbor.hpp"

std::vector<std::vector<int>> neighbor::direct_find(const int np, const std::vector<double> &sp,
                                                    const std::vector<double> &xp, const std::vector<double> &yp, const int neighbor_scale)
{
    /*************************************************************************
    **    Subroutine to calculate the smoothing funciton for each particle and
    **    the interaction parameters used by the SPH algorithm. Interaction 
    **    pairs are determined by directly comparing the particle distance 
    **    with the corresponding smoothing length.
    **      itimestep : Current time step                                 [in]
    **      np        : Number of particles                               [in]
    **      sp        : Smoothing Length                                  [in]
    **      x         : Coordinates of all particles                      [in]
    **      niac      : Number of interaction pairs                      [out]
    **      pair_i    : List of first partner of interaction pair        [out]
    **      pair_j    : List of second partner of interaction pair       [out]
    **      w         : Kernel for all interaction pairs                 [out]
    **      dwdx      : Derivative of kernel with respect to x, y and z  [out]
    **      countiac  : Number of neighboring particles                  [out]
    *************************************************************************/
    // ** internal variables
    // int sumiac, maxiac, miniac, noiac, maxp, minp;
    double _dx, _dy, _dr, _msp;
    std::vector<std::vector<int>> _neighborlist(np, std::vector<int>()); // ! returning values
    // int *countiac = new int [np];
    //// if      (skf == 1){neighbor_scale = 2;}
    //// else if (skf == 2){neighbor_scale = 3;}/* default: neighbor_scale = 3 */
    //// else if (skf == 3){neighbor_scale = 3;}

    ////pair_i.clear();
    ////pair_j.clear();
    for (int i = 0; i < np - 1; i++)
    {
        for (int j = i + 1; j < np; j++)
        {
            _dx = xp[i] - xp[j];
            _dy = yp[i] - yp[j];
            _dr = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2));
            _msp = (sp[i] + sp[j]) / 2;
            // msp = sp[i] < sp[j] ? sp[i] : sp[j]; // ! for DC PSE's neighbour searching
            if (_dr < neighbor_scale * _msp)
            {
                // // pair_i.push_back(i);
                // // pair_j.push_back(j);
                _neighborlist[i].push_back(j);
                _neighborlist[j].push_back(i);
                // countiac[i] += 1;
                // countiac[j] += 1;
            }
        }
    }

    // // statisticas for the interaction
    // sumiac = 0;
    // maxiac = 0;
    // miniac = 1000;
    // noiac  = 0;
    // for (int i = 0; i < np; i++)
    // {
    //     sumiac = sumiac + countiac[i];
    //     if (countiac[i] > maxiac)
    //     {
    //         maxiac    = countiac[i];
    //         maxp      = i;
    //     }
    //     else if (countiac[i] < miniac)
    //     {
    //         miniac    = countiac[i];
    //         minp      = i;
    //     }
    //     else if (countiac[i] == 0)
    //     {
    //         noiac = noiac + 1;
    //     }
    // }

    // // -- delete internal variables
    // delete[] countiac;

    //// printf("\n >> Statistics: interactions per particle:");
    //// printf("\n**** Particle: %d maximal interactions: %d", maxp, maxiac);
    //// printf("\n**** Particle: %d minimal interactions: %d", minp, miniac);
    //// printf("\n**** Average : %f", float(sumiac) / float(np) );
    //// printf("\n**** Total pairs : %d", niac);
    //// printf("\n**** Particles with no interactions: %d \n", noiac);

    return _neighborlist;
}

void neighbor::direct_3d(Particle &p, const int np, const double neighbor_scale){
    for (int i = 0; i < np - 1; i++)
    {
        for (int j = i + 1; j < np; j++)
        {
            double _dx = p.x[i] - p.x[j];
            double _dy = p.y[i] - p.y[j];
            double _dz = p.z[i] - p.z[j];
            double _dr = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2) + std::pow(_dz, 2));
            double _msp = (p.s[i] + p.s[j]) / 2;

            if (_dr < neighbor_scale)
            {
                p.neighbor[i].push_back(j);
                p.neighbor[j].push_back(i);
            }
        }
    }
}
