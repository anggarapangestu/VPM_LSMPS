#include "neighbor.hpp"

// std::vector<std::vector<int>> hash(const int np, const std::vector<double> &sp,
// 										   const std::vector<double> &xp, const std::vector<double> &yp, int neighbor_scale)
// {
//     for (size_t i = 0; i < np; i++){

//     }
// }

void search_cell(const int i, Particle &p, Cell &cell, int loc, int loc_x, int loc_y, int loc_z, double range){
    
    int cell_num = (loc_x) + (loc_y) * 4 + (loc_z) * 4 * 4;

    for(int j = 0; j < cell.particle_inside[cell_num].size(); j++){

        int _n = cell.particle_inside[cell_num][j];
        
        double r = std::sqrt(std::pow(p.x[loc] - p.x[_n], 2) + std::pow(p.y[loc] - p.y[_n] , 2) + std::pow(p.z[loc] - p.z[_n] , 2));

        if (r < range && r != 0)
            p.neighbor[i].push_back(_n);
    }
}

void neighbor::hash(Particle &p, Cell &c, const int np, double cell_size)
{
    clock_t t = clock();
    int _x[27] = {0, 1, -1, 0,  0, 0,  0, 0,  0,  0,  0, 1, -1,  1, -1, 1,  1, -1, -1, 1,  1,  1,  1, -1, -1, -1, -1};
    int _y[27] = {0, 0,  0, 1, -1, 0,  0, 1,  1, -1, -1, 1,  1, -1, -1, 0,  0,  0,  0, 1,  1, -1, -1,  1,  1, -1, -1};
    int _z[27] = {0, 0,  0, 0,  0, 1, -1, 1, -1,  1, -1, 0,  0,  0,  0, 1, -1,  1, -1, 1, -1,  1, -1,  1, -1, -1,  1};
    
    int cell_width = (Pars::xmax - Pars::xmin)/cell_size; //Saat ini domain masih berbentuk kubus
    
    for (int i = 0; i < np; i++){
        int loc = p.hash_cell[i];
        //printf("<%d, %d, %d> \n", c.x[loc], c.y[loc], c.z[loc]);
        for (int j = 0; j < 27; j++){
            if(((c.x[loc] + _x[j] <= (cell_width - 1)) && (c.x[loc] + _x[j] >= 0)) && 
               ((c.y[loc] + _y[j] <= (cell_width - 1)) && (c.y[loc] + _y[j] >= 0)) && 
               ((c.z[loc] + _z[j] <= (cell_width - 1)) && (c.z[loc] + _z[j] >= 0))){
                    //3.5* particle size is effective radius for LSMPS (see: LSMPS 2014)  
                    search_cell(i, p, c, i, c.x[loc] + _x[j], c.y[loc] + _y[j], c.z[loc] + _z[j], 3.5 * p.s[i]);
                    //printf("oke\n");
            }
        }
    }
    t = clock() - t;
	printf("Neighbor search using Spatial Hash [%fs]\n", (double)t/CLOCKS_PER_SEC);
}