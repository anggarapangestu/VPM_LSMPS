#include "neighbor.hpp"

void neighbor::inter_search(const Particle &p_src, const Particle &p_tar, std::vector<std::vector<int>>& ngh_list)
{
	// Initialize the target and source position 
	// *****************************************
	// Target for neighbor list
	const std::vector<double> &x_tar = p_tar.x;
	const std::vector<double> &y_tar = p_tar.y;
	const std::vector<double> &s_tar = p_tar.s;
	// Source of particle list
	const std::vector<double> &x_src = p_src.x;
	const std::vector<double> &y_src = p_src.y;
	const std::vector<double> &s_src = p_src.s;

	// const static int R_e = 3.10e0; // ! should be change as class internal variables

	// Create a hash box container parameter
	// *************************************
	double *_mindom = new double[Pars::DIM]; // minimum particle location
	double *_maxdom = new double[Pars::DIM]; // maximum particle location

	double _sp_min = *std::min_element(s_src.begin(), s_src.end()); // minimum core size of scattered particle
	double _sp_max = *std::max_element(s_src.begin(), s_src.end()); // maximum core size of scattered particle
	double _sp_av = (_sp_min + _sp_max) / 2;						// average core size of scattered particle

	_mindom[0] = *std::min_element(x_src.begin(), x_src.end()) - Pars::r_sup * _sp_av;
	_mindom[1] = *std::min_element(y_src.begin(), y_src.end()) - Pars::r_sup * _sp_av;
	_maxdom[0] = *std::max_element(x_src.begin(), x_src.end()) + Pars::r_sup * _sp_av;
	_maxdom[1] = *std::max_element(y_src.begin(), y_src.end()) + Pars::r_sup * _sp_av;

	// Create the hash box containing the source particle
	double _dbox = Pars::r_sup * Pars::r_buff * _sp_max; 	// size of single box
	int *_nbox = new int[2];			   	// store number of boxes
	for (size_t i = 0; i < Pars::DIM; i++)
	{
		_nbox[i] = std::ceil((_maxdom[i] - _mindom[i]) / _dbox);
	}

	// Create the hash bucket
	// **********************
	std::vector<std::vector<std::vector<int>>> _neighborParBox(_nbox[0], std::vector<std::vector<int>>(_nbox[1], std::vector<int>()));
	// Put all the source particle into the occupied hash bucket 
	for (size_t i = 0; i < x_src.size(); i++)
	{
		int _idx_i = std::floor((x_src[i] - _mindom[0]) / _dbox);
		int _idx_j = std::floor((y_src[i] - _mindom[1]) / _dbox);
		_neighborParBox[_idx_i][_idx_j].push_back(i);
	}

	// Evaluate the neighbor
	// *********************
	// Resize the neighbor list array
	ngh_list.clear();
	ngh_list.resize(x_tar.size(), std::vector<int>());
	// Evaluate to all particle in the box
	for (size_t i = 0; i < x_tar.size(); i++)	// Iterate through each target particle
	{
		// Determine the box position of the current target particle located 
		int _idx_i = std::floor((x_tar[i] - _mindom[0]) / _dbox);
		int _idx_j = std::floor((y_tar[i] - _mindom[1]) / _dbox);

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

					// Neighbor evaluation
					_dx = x_src[pot_ngh_ID] - x_tar[i];
					_dy = y_src[pot_ngh_ID] - y_tar[i];
					_dr = std::sqrt(_dx*_dx + _dy*_dy);
					// _sij = (s_tar[i] + s_src[pot_ngh_ID]) / 2.0e0;
					_sij = s_tar[i];
					if (_dr <= Pars::r_sup * Pars::r_buff * _sij)
					{
						ngh_list[i].push_back(pot_ngh_ID);
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