#include "neighbor.hpp"
#include "base_grid.hpp"

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
    double *dx = new double[Pars::DIM];
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
    delete[] countiac;

    return _neighborlist;
}


void base_grid::create_grid(double xfixed, double hG, double xGmin, double xGmax, 
  int nG, int &nG_new, std::vector<double> &xGi)
{
	// Internal variables 
	int nG1, nG2;

  if ( (xGmin <= xfixed) && (xfixed <= xGmax) )
  {
    nG1 = std::ceil( (xfixed - xGmin)/hG); // Number of nodes from fixed point to the min point
    nG2 = std::ceil( (xGmax - xfixed)/hG); // Number of nodes from fixed point to the max point
    nG_new = nG1 + nG2 + 1;
    xGi.resize(nG_new);

    for (int i = 0; i < nG1; i++){
      // printf("\n part one %d %d", i, nG1);
      xGi[i] = xfixed - (nG1 - i)*hG;
    }
    // printf("done");
    for (int i = 0; i <= nG2; i++){
      // printf("\n%d %d", i, nG2);
      xGi[nG1+i] = xfixed + i*hG;
    }
  }
  else if (xfixed < xGmin)
  {
    nG1 = std::floor( (xGmin - xfixed)/hG); // Number of nodes from fixed point to the min point
    nG2 = std::ceil ( (xGmax - xfixed)/hG); // Number of nodes from fixed point to the max point
    for (int i = 0; i < (nG2-nG1); i++)
      xGi[i] = xfixed + (nG1 + i)*hG;
    nG_new = (nG2-nG1);
    xGi.resize(nG_new);
    // nG_new = nG; 
  }
  else if (xGmax < xfixed)
  {
    nG1 = std::ceil ( (xfixed - xGmin)/hG); // Number of nodes from fixed point to the min point
    nG2 = std::floor( (xfixed - xGmax)/hG); // Number of nodes from fixed point to the max point
    for (int i = 0; i < (nG1-nG2); i++)
      xGi[i] = xfixed - (nG2 + i)*hG;
    nG_new = (nG1-nG2);
    xGi.resize(nG_new);
    // nG_new = nG;
  }
}

void base_grid::find_cell(double xPi, double yPi, int nGx, int nGy, 
    const std::vector<double> &xGi, const std::vector<double> &yGi, double hG, int &i_cell, int &j_cell)
{
    // ===================================================================
    // Finding the cell (named by it's left-side) containing particle i
    //                         .pi
    //                                  .pi
    //                          |__________|
    //                          |          |         
    //                          |     .pi  |
    //                          |__________|
    //                          |          |
    //                        left     left+1=right
    // ===================================================================
    // internal variables
    double xminG, yminG;
    // find the left side:
    //    x-dir
    xminG   = *std::min_element(xGi.begin(), xGi.begin() + nGx);
    i_cell  = std::floor((xPi - xminG) / hG);
    //    y-dir
    yminG   = *std::min_element(yGi.begin(), yGi.begin() + nGy);
    j_cell  = std::floor((yPi - yminG) / hG);
}

void base_grid::find_neighborcells(int i_cell, int j_cell, int neighbor_scale, int nGx, int nGy,
    int &lx_nbc, int &rx_nbc, int &ly_nbc, int &ry_nbc)
{
      // ---define the west (x_loc -2) and east (x_loc+3)
    lx_nbc = i_cell - neighbor_scale;
    rx_nbc = i_cell + neighbor_scale;
    // ---define the south (y_loc -2) and north (y_loc+3)
    ly_nbc = j_cell - neighbor_scale;
    ry_nbc = j_cell + neighbor_scale;
    
    // // !!!===================================================================
    // // ! Finding the neighbor cells around the particle -i
    // // !           ____ _____________________________
    // // !          |    |    |    |    |    |    |    |
    // // !          |____|____|____|____|____|____|____|
    // // !          |    |    |    |    |    |    |    |
    // // !          |____|____|____|____|____|____|____|
    // // !          |    |    |    |    |    |    |    |
    // // !          |____|____|____|____|____|____|____|
    // // !          |    |    |    | .p |    |    | p2.|
    // // !          |____|____|____|____|____|____|____|
    // // !          |    |    |    |    |    |    |    |
    // // !          |____|____|____|____|____|____|____|
    // // !          |    |    |    |    |    |    |    |
    // // !          |____|____|____|____|____|____|____|
    // // !          |    |    |    |    |    |    |    |
    // // !          |____|____|____|____|____|____|____|
    // // !left_nbc=-3   -2   -1  left   +1  +2   +3   +4 = right_nbc
    // // ! We need information input: from find_cell  we have name of cell containning the particle i_cell,j_cell,k_cell
    // // ! We need to know how many neighbor cells around, we want to take, controled by: neighbor_scale --> WE CAN CHOSE BIGER THAN WHAT WE NEED 
    // // !           nb_scale is |x| or |xpi-xGi|/hG tell that distant between particle and node/cell = how many "time" of Grid spacing
    // // ! NOTE THAT: - for finding nodes arround particle in a cell we just need  (i_cell-nb_scale) --> (i_cell+nb_scale+1). e.g. for P2g, G2P
    // // !            - for finding "particles in cell" around we just need  (i_cell-nb_scale) --> (i_cell+nb_scale)  . i.e. for link_list
    // // !            - BUT for easier we let all with (i_cell-nb_scale) --> (i_cell+nb_scale), actually add more 1 or not, not really affect to the result:
    // // ! e.x. nb_scale=3; if p close to left, then may "-3" affect to p, if p close to "+1", then may "+4" affects to p. That is MAYBE
    // // ! bcz that happend when p lies on "+1" the the cell now not is "left" any more, that is "+1" , and +4 =+1+3. Assume that happen then kernal function at that node ~> 0
    // // ! But for searching particle , may the p2 (cell of +3) close to +4 can affect to p
    // // ! Finaly we could obtain: left and right neighbor cell (l_nbc, r_nbc)
    // // !!!===================================================================

    // // !!!=== internal ===
     
    // // !=== Finding the nieghboring nodes arround the particle (i,j) (moved)
    // // ! x view
    // if (i_cell <= neighbor_scale)
    // {
    //     lx_nbc = 0;
    //     rx_nbc = i_cell + neighbor_scale;
    // }
    // else if (i_cell >= (nGx-neighbor_scale))
    // {
    //     lx_nbc = i_cell-neighbor_scale;
    //     rx_nbc = nGx;
    // }
    // else
    // {
    //     lx_nbc = i_cell-neighbor_scale;
    //     rx_nbc = i_cell+neighbor_scale;
    // }
    // // ! y view
    // if (j_cell <= neighbor_scale)
    // {
    //     ly_nbc = 0;
    //     ry_nbc = j_cell+neighbor_scale;
    // }
    // else if (j_cell >= (nGy-neighbor_scale))
    // {
    //     ly_nbc = j_cell-neighbor_scale;
    //     ry_nbc = nGy;
    // }
    // else
    // {
    //     ly_nbc = j_cell-neighbor_scale;
    //     ry_nbc = j_cell+neighbor_scale;
    // }
}