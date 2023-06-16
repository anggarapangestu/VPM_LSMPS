#include "base_grid.hpp"

void hashing(const int np, const std::vector<double> &xp, const std::vector<double> &yp, 
             const std::vector<double> &zp, std::vector<int> &hx, std::vector<int> &hy,
             std::vector<int> &hz, std::vector<std::vector<int>> &hinside, double cell_size)
{
    //initialize values;
    hx.clear(); hx.resize(np);
    hy.clear(); hy.resize(np);
    hz.clear(); hz.resize(np);
    hinside.clear(); hinside.resize(np);
    
    double conversionfactor = 1.0 / cell_size; //stil not fixed yet
    
    for (int i = 0 ; i < np; i++){
        int cell_num = floor(xp[i] * conversionfactor) + floor(yp[i] * conversionfactor) * cell_size + floor(zp[i] * conversionfactor) * cell_size * cell_size;
        hinside[cell_num].push_back(i);
        hx.push_back(floor(xp[i] * conversionfactor));
        hy.push_back(floor(yp[i] * conversionfactor));
        hz.push_back(floor(zp[i] * conversionfactor));

    }
}

void base_grid::spatial_hashing(Particle &p, Cell &cell, const int np, const double cell_size){

    cell.x.clear(); cell.x.resize(np);
    cell.y.clear(); cell.y.resize(np);
    cell.z.clear(); cell.z.resize(np);
    cell.particle_inside.clear(); cell.particle_inside.resize(4*4*4);
    p.hash_cell.clear();

    //Particle _par;
    //_par = p; 
    double xmin = *std::min_element(p.x.begin(), p.x.end());
    double ymin = *std::min_element(p.y.begin(), p.y.end());
    double zmin = *std::min_element(p.z.begin(), p.z.end());
    double xmax = *std::max_element(p.x.begin(), p.x.end());

    double cell_width = (xmax - xmin)/cell_size; //saat ini domain masih berbentuk kubus
    double conversionfactor = 1.0 / cell_size;
    
    
    for (int i = 0 ; i < np; i++){        
        int cell_num = floor((p.x[i] + abs(xmin)) * conversionfactor) + floor((p.y[i] + abs(ymin)) * conversionfactor) * cell_width + floor((p.z[i] + abs(zmin)) * conversionfactor) * cell_width * cell_width;
        cell.particle_inside[cell_num].push_back(i);
        cell.x[i] = (floor((p.x[i] + abs(xmin)) * conversionfactor));
        cell.y[i] = (floor((p.y[i] + abs(ymin)) * conversionfactor));
        cell.z[i] = (floor((p.z[i] + abs(ymin)) * conversionfactor));
        p.hash_cell.push_back(cell_num);
    }
}