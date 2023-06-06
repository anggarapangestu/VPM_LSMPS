#include "neighbor.hpp"

std::vector<std::vector<int>> neighbor::link_list(const int np, const std::vector<double> &sp,
                                                  const std::vector<double> &xp, const std::vector<double> &yp, int neighbor_scale)
{
    // ----------------------------------------------------------------------
    //    Subroutine to calculate the smoothing funciton for each particle and
    //    the interaction parameters used by the SPH algorithm. Interaction
    //    pairs are determined by using a sorting grid linked list
    //    the complexity of the linked-list algorithm is of order O(N).
    //    (Computation time is of order O(N))
    //  np             : Number of particles                                 [in]
    //  sp             : Smoothing Length/core-size/grid-spacing             [in]
    //  x              : Coordinates of all particles                        [in]
    //  neighbor_scale : how many time of Smoothing Length want to take      [in]
    //  pair_i         : List of first partner of interaction pair          [out]
    //  pair_j         : List of second partner of interaction pair         [out]
    //  countiac       : Number of neighboring particles                    [out]

    // Internal variables
    double grid_space;
    int i, j, d, ix, jx;
    int xcell, ycell;
    double dr, msp;
    int *celldata = new int[np];
    int *minxcell = new int[2];
    int *maxxcell = new int[2];
    double *dx = new double[Pars::dim];
    int *i_c = new int[np];
    int *j_c = new int[np];
    std::vector<std::vector<int>> _neighborlist(np, std::vector<int>()); // ! returning values

    // For Statistics
    int ngx_ll, ngy_ll;
    int pr, sumiac, maxiac, noiac, miniac, maxp, minp;
    int *countiac = new int[np];
    double *xgmax = new double[2];
    double *xgmin = new double[2];
    int *ncgrid = new int[2];

    // Counting number of interaction of each particle
    // Initialize grid:
    //    Mesh spacing: Maximum Smoothing length of all particles/ max_core-size: take maximum to satisfy taking neighbor_scale*mesh_spacing around,
    //    if mesh_spacing small that lead to missing particles arounding, bcz neighbor_scale*mesh_spacing becomes too small for bigger sp
    grid_space = *std::max_element(sp.begin(), sp.begin() + np); // grid spacing
    //    Range of sorting grid
    xgmax[0] = *std::max_element(xp.begin(), xp.end()) + (neighbor_scale + 1) * grid_space; //    add 2 more node (each node for each end side)
    xgmax[1] = *std::max_element(yp.begin(), yp.end()) + (neighbor_scale + 1) * grid_space; //    add 2 more node (each node for each end side)
    xgmin[0] = *std::min_element(xp.begin(), xp.end()) - (neighbor_scale + 1) * grid_space;
    xgmin[1] = *std::min_element(yp.begin(), yp.end()) - (neighbor_scale + 1) * grid_space;
    
    // Calculate grid size first for saving initial generating memory:
    ncgrid[0] = std::ceil((xgmax[0] - xgmin[0]) / grid_space) + 1; // add 2 more node
    ncgrid[1] = std::ceil((xgmax[1] - xgmin[1]) / grid_space) + 1;
    
    // to make sure that: there is no particle lies on sides of box, and all particles inside the box
    // ifelse Error of Sofware Fortran lead to wrong result
    // === allocate redistribution grids ===
    std::vector<double> xgrid(ncgrid[0] + 2 + 3);
    std::vector<double> ygrid(ncgrid[1] + 2 + 3);
    d_base_grid.create_grid(0.e0, grid_space, xgmin[0], xgmax[0], ncgrid[0], ngx_ll, xgrid);
    d_base_grid.create_grid(0.e0, grid_space, xgmin[1], xgmax[1], ncgrid[1], ngy_ll, ygrid);

    //    Initialize grid ---> HEAD of Chain
    std::vector<std::vector<double>> grid(ngx_ll, std::vector<double>(ngy_ll));
    
    //    Position particles on grid and create linked list:
    for (int i = 0; i < np; i++)
    {
        d_base_grid.find_cell(xp[i], yp[i], ngx_ll, ngy_ll, xgrid, ygrid, grid_space, i_c[i], j_c[i]);
        celldata[i] = grid[i_c[i]][j_c[i]]; // --> head chain =0 --> 1st i --> 2nd i (same cell) ...link list.. --> before "end of chain"
        grid[i_c[i]][j_c[i]] = i;           // --> end of chain
    }

    xgrid.clear();
    ygrid.clear();
    //    Determine interaction parameters:
    //// pair_i.clear();
    //// pair_j.clear();
    for (int i = 0; i < np; i++)
    {
        //      Determine range of grid to go through:
        d_base_grid.find_neighborcells(i_c[i], j_c[i], neighbor_scale, ngx_ll, ngy_ll, minxcell[0], maxxcell[0], minxcell[1], maxxcell[1]);
        //      Search grid:
        for (int ycell = minxcell[1]; ycell <= maxxcell[1]; ycell++)
        {
            for (int xcell = minxcell[0]; xcell <= maxxcell[0]; xcell++)
            {
                j = grid[xcell][ycell];
                while (j > i)
                {
                    double _dx = xp[i] - xp[j];
                    double _dy = yp[i] - yp[j];
                    double _dr = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2));
                    // -----Symmetrization of particle interaction:
                    msp = sp[i]; // (sp[i] + sp[j]) / 2.0e0;
                    if (_dr < (((double)neighbor_scale + 0.1) * msp))
                    {
                        // Neighboring pair list, and total interaction number and
                        // the interaction number for each particle
                        //// pair_i.push_back(i);
                        //// pair_j.push_back(j);
                        _neighborlist[i].push_back(j);
                        _neighborlist[j].push_back(i);
                        // countiac[i] += 1;
                        // countiac[j] += 1;
                    }
                    j = celldata[j]; //   down one level to the previous particle (smaller) which is in cell i as well
                }
            }
        }
    }
    // -- delete internal variables
    for (size_t i = 0; i < ngx_ll; i++)
    {
        grid[i].clear();
    }
    grid.clear();

    delete[] celldata;
    delete[] minxcell;
    delete[] maxxcell;
    delete[] dx;
    delete[] i_c;
    delete[] j_c;
    delete[] xgmax;
    delete[] xgmin;
    delete[] ncgrid;

    // // Statistics for the interaction
    // sumiac = 0;
    // maxiac = 0;
    // miniac = 1000;
    // noiac = 0;
    // for (int i = 0; i < np; i++)
    // {
    //     sumiac = sumiac + countiac[i];
    //     if (countiac[i] > maxiac)
    //     {
    //         maxiac = countiac[i];
    //         maxp = i;
    //     }
    //     else if (countiac[i] < miniac)
    //     {
    //         miniac = countiac[i];
    //         minp = i;
    //     }
    //     else if (countiac[i] == 0)
    //     {
    //         noiac = noiac + 1;
    //     }
    // }
    delete[] countiac;

    //// printf("\n >> Statistics: interactions per particle:");
    //// printf("\n**** Particle: %d maximal interactions: %d", maxp, maxiac);
    //// printf("\n**** Particle: %d minimal interactions: %d", minp, miniac);
    //// printf("\n**** Average : %f", float(sumiac) / float(np) );
    //// printf("\n**** Total pairs : %d", niac);
    //// printf("\n**** Particles with no interactions: %d \n", noiac);

    return _neighborlist;
}
