#include "neighbor.hpp"

void neighbor::spatial_hash(Particle &parEval){
    // Initialize the position 
	// *****************************************
	const std::vector<double> &_x = parEval.x;
	const std::vector<double> &_y = parEval.y;
	const std::vector<double> &_s = parEval.s;

	// const static int R_e = 3.10e0; // ! should be change as class internal variables

	// Create a hash box container parameter
	// *************************************
	double *_mindom = new double[Pars::DIM]; // minimum particle location
	double *_maxdom = new double[Pars::DIM]; // maximum particle location

	double _sp_min = *std::min_element(_s.begin(), _s.end()); // minimum core size of scattered particle
	double _sp_max = *std::max_element(_s.begin(), _s.end()); // maximum core size of scattered particle
	double _sp_av = (_sp_min + _sp_max) / 2;						// average core size of scattered particle

	_mindom[0] = *std::min_element(_x.begin(), _x.end()) - Pars::r_sup * _sp_av;
	_mindom[1] = *std::min_element(_y.begin(), _y.end()) - Pars::r_sup * _sp_av;
	_maxdom[0] = *std::max_element(_x.begin(), _x.end()) + Pars::r_sup * _sp_av;
	_maxdom[1] = *std::max_element(_y.begin(), _y.end()) + Pars::r_sup * _sp_av;

	// Create the hash box containing the source particle
	double _dbox = Pars::r_sup * Pars::r_buff * _sp_av; // size of single box
	int *_nbox = new int[2];			   // store number of boxes
	for (size_t i = 0; i < Pars::DIM; i++)
	{
		_nbox[i] = std::ceil((_maxdom[i] - _mindom[i]) / _dbox);
	}

	// Create the hash bucket
	// **********************
	std::vector<std::vector<std::vector<int>>> _neighborParBox(_nbox[0], std::vector<std::vector<int>>(_nbox[1], std::vector<int>()));
	// Put all the source particle into the occupied hash bucket 
	for (size_t i = 0; i < _x.size(); i++)
	{
		int _idx_i = std::floor((_x[i] - _mindom[0]) / _dbox);
		int _idx_j = std::floor((_y[i] - _mindom[1]) / _dbox);
		_neighborParBox[_idx_i][_idx_j].push_back(i);
	}

	// Evaluate the neighbor
	// *********************
	// Resize the neighbor list array
	parEval.neighbor.clear();
	parEval.neighbor.resize(parEval.num, std::vector<int>());
	// Evaluate to all particle in the box
	for (size_t i = 0; i < parEval.num; i++)	// Iterate through each target particle
	{
		// Determine the box position of the current target particle located 
		int _idx_i = std::floor((_x[i] - _mindom[0]) / _dbox);
		int _idx_j = std::floor((_y[i] - _mindom[1]) / _dbox);

		// Internal variable
		// std::vector<int> _parBox{};
		int ix_start, ix_end;
		int iy_start, iy_end;
		int pot_ngh_ID;
		double _dx, _dy, _dr, _sij;

		// Set the location of neighboring box
		ix_start = _idx_i - 1;
		ix_end = _idx_i + 1;
		iy_start = _idx_j - 1;
		iy_end = _idx_j + 1;

		// Evaluate each variable 
		if (ix_start < 0){
			ix_start = 0;
		}
		if (ix_end > _nbox[0] - 1){
			ix_end = _nbox[0] - 1;
		}
		if (iy_start < 0){
			iy_start = 0;
		}
		if (iy_end < _nbox[1] - 1){
			iy_end = _nbox[1] - 1;
		}

		// Evaluate the neighbor with the particle in each neighboring box
		for (size_t ix = ix_start; ix <= ix_end; ix++)
		{
			for (size_t iy = iy_start; iy <= iy_end; iy++)
			{
				for (size_t iPar = 0; iPar < _neighborParBox[ix][iy].size(); iPar++)
				{
					// Take the current potential neighbor ID
					pot_ngh_ID = _neighborParBox[ix][iy][iPar];
                    
                    // No repeating ID
                    if (pot_ngh_ID == i)
                    {
                        continue;
                    }

					// Neighbor evaluation
					_dx = _x[pot_ngh_ID] - _x[i];
					_dy = _y[pot_ngh_ID] - _y[i];
					_dr = std::sqrt(_dx*_dx + _dy*_dy);
					// _sij = (_s[i] + _s[pot_ngh_ID]) / 2.0e0;
					_sij = _s[i];
					if (_dr <= Pars::r_sup * Pars::r_buff * _sij)
					{
						parEval.neighbor[i].push_back(pot_ngh_ID);
					}

				}
			}
		}
	}
	
	// Free the space
	delete[] _mindom;
	delete[] _maxdom;
	delete[] _nbox;

	return;
}

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
    
    double xmin = *std::min_element(p.x.begin(), p.x.end());
	double xmax = *std::max_element(p.x.begin(), p.x.end());
	int cell_width = (xmax - xmin)/cell_size; //Saat ini domain masih berbentuk kubus
    
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