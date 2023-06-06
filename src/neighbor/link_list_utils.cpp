#include "base_grid.hpp"

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
