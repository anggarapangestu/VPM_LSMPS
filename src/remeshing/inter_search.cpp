#include "base_remeshing.hpp"

std::vector<std::vector<int>> base_remeshing::inter_search(const Particle &p_scatter, const Particle &p_grid)
{
	// -- generate box
	const std::vector<double> &_x_grd = p_grid.x;
	const std::vector<double> &_y_grd = p_grid.y;
	const std::vector<double> &_s_grd = p_grid.s;
	const std::vector<double> &_x_par = p_scatter.x;
	const std::vector<double> &_y_par = p_scatter.y;
	const std::vector<double> &_s_par = p_scatter.s;

	const static int R_e = 3.10e0; // ! should be change as class internal variables

#pragma region multibox
	double *_mindom = new double[Pars::dim]; // minimum particle location
	double *_maxdom = new double[Pars::dim]; // maximum particle location

	double _sp_min = *std::min_element(_s_par.begin(), _s_par.end()); // minimum core size of scattered particle
	double _sp_max = *std::max_element(_s_par.begin(), _s_par.end()); // maximum core size of scattered particle
	double _sp_av = (_sp_min + _sp_max) / 2;						  // average core size of scattered particle

	_mindom[0] = *std::min_element(_x_par.begin(), _x_par.end()) - Pars::r_scale * _sp_av;
	_mindom[1] = *std::min_element(_y_par.begin(), _y_par.end()) - Pars::r_scale * _sp_av;
	_maxdom[0] = *std::max_element(_x_par.begin(), _x_par.end()) + Pars::r_scale * _sp_av;
	_maxdom[1] = *std::max_element(_y_par.begin(), _y_par.end()) + Pars::r_scale * _sp_av;

	double _dbox = Pars::r_scale * _sp_av; // size of single box
	int *_nbox = new int[2];			   // store number of boxes
	for (size_t i = 0; i < Pars::dim; i++)
	{
		_nbox[i] = std::ceil(std::abs(_maxdom[i] - _mindom[i]) / _dbox);
	}

	// -- assign neighborhood
	std::vector<std::vector<std::vector<int>>> _neighborParBox(_nbox[0], std::vector<std::vector<int>>(_nbox[1], std::vector<int>()));
	// _neighborParBox: The spatial table to store all particle index in the correspond grid
	for (size_t i = 0; i < _x_par.size(); i++)
	{
		int _idx_i = std::floor((_x_par[i] - _mindom[0]) / _dbox);
		int _idx_j = std::floor((_y_par[i] - _mindom[1]) / _dbox);
		_neighborParBox[_idx_i][_idx_j].push_back(i);
	}

	// -- search neighbor of the grid
	std::vector<std::vector<int>> _neighborGrid(_x_grd.size(), std::vector<int>());
	for (size_t i = 0; i < _x_grd.size(); i++)
	{
		if (_x_grd[i] >= _mindom[0] && _x_grd[i] <= _maxdom[0] &&
			_y_grd[i] >= _mindom[1] && _y_grd[i] <= _maxdom[1])
		{
			int _idx_i = std::floor((_x_grd[i] - _mindom[0]) / _dbox);
			int _idx_j = std::floor((_y_grd[i] - _mindom[1]) / _dbox);

			std::vector<int> _parBox{};
			int ix_start, ix_end;
			int iy_start, iy_end;
	
			if (_idx_i > 0 && _idx_i < (_nbox[0] - 1) &&
				_idx_j > 0 && _idx_j < (_nbox[1] - 1))
			{
				ix_start = _idx_i - 1;
				ix_end = _idx_i + 1;
				iy_start = _idx_j - 1;
				iy_end = _idx_j + 1;
			}
			else if (_idx_i == 0 &&
					 _idx_j > 0 && _idx_j < (_nbox[1] - 1))
			{
				ix_start = _idx_i;
				ix_end = _idx_i + 1;
				iy_start = _idx_j - 1;
				iy_end = _idx_j + 1;
			}
			else if (_idx_i == (_nbox[0] - 1) &&
					 _idx_j > 0 && _idx_j < (_nbox[1] - 1))
			{
				ix_start = _idx_i - 1;
				ix_end = _idx_i;
				iy_start = _idx_j - 1;
				iy_end = _idx_j + 1;
			}
			else if (_idx_j == 0 &&
					 _idx_i > 0 && _idx_i < (_nbox[0] - 1))
			{
				ix_start = _idx_i - 1;
				ix_end = _idx_i + 1;
				iy_start = _idx_j;
				iy_end = _idx_j + 1;
			}
			else if (_idx_j == (_nbox[1] - 1) &&
					 _idx_i > 0 && _idx_i < (_nbox[0] - 1))
			{
				ix_start = _idx_i - 1;
				ix_end = _idx_i + 1;
				iy_start = _idx_j - 1;
				iy_end = _idx_j;
			}
			else if (_idx_i == 0 && _idx_j == 0)
			{
				ix_start = _idx_i;
				ix_end = _idx_i + 1;
				iy_start = _idx_j;
				iy_end = _idx_j + 1;
			}
			else if (_idx_i == (_nbox[0] - 1) && _idx_j == (_nbox[1] - 1))
			{
				ix_start = _idx_i - 1;
				ix_end = _idx_i;
				iy_start = _idx_j - 1;
				iy_end = _idx_j;
			}
			else if (_idx_i == (_nbox[0] - 1) && _idx_j == 0)
			{
				ix_start = _idx_i - 1;
				ix_end = _idx_i;
				iy_start = _idx_j;
				iy_end = _idx_j + 1;
			}
			else if (_idx_i == 0 && _idx_j == (_nbox[1] - 1))
			{
				ix_start = _idx_i;
				ix_end = _idx_i + 1;
				iy_start = _idx_j - 1;
				iy_end = _idx_j;
			}
			// Input tetangga dari box ke i 
			for (size_t ix = ix_start; ix <= ix_end; ix++)
			{
				for (size_t iy = iy_start; iy <= iy_end; iy++)
				{
					for (size_t iPar = 0; iPar < _neighborParBox[ix][iy].size(); iPar++)
					{
						_parBox.push_back(_neighborParBox[ix][iy][iPar]);
					}
				}
			}

			// std::vector<int> _parBox = _neighborParBox[_idx_i][_idx_j];
			for (size_t j = 0; j < _parBox.size(); j++)
			{
				double _dx = _x_par[_parBox[j]] - _x_grd[i];
				double _dy = _y_par[_parBox[j]] - _y_grd[i];
				double _dr = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2));
				double _sij = (_s_grd[i] + _s_par[_parBox[j]]) / 2.0e0;
				if (_dr <= R_e * _sij)
				{
					_neighborGrid[i].push_back(_parBox[j]);
				}
			}
		}
		else
		{
			continue;
		}
	}

	delete[] _mindom;
	delete[] _maxdom;
	delete[] _nbox;
#pragma endregion

	// std::vector<std::vector<int>> _neighborGridTemp(_x_grd.size(), std::vector<int>());
	// #pragma omp parallel for
	// for (int i = 0; i < _x_grd.size(); i++)
	// {
	// 	for (size_t j = 0; j < _x_par.size(); j++)
	// 	{
	// 		double _dr = std::sqrt(std::pow((_x_par[j] - _x_grd[i]), 2) + std::pow((_y_par[j] - _y_grd[i]), 2));
	// 		double _sp = (_s_grd[i] + _s_par[j]) * 0.5;
	// 		if (_dr <= R_e * _sp)
	// 		{
	// 			// if (i != j)
	// 			// {
	// 				_neighborGridTemp[i].push_back(j);
	// 			// }
	// 		}
	// 	}
	// }

	// for (size_t i = 0; i < _x_grd.size(); i++)
	// {
	// 	if (!_neighborGrid[i].empty())
	// 	{
	// 		std::cout << i << " " << _neighborGrid[i].size() << " " << _neighborGridTemp[i].size() << std::endl;
	// 	}
	// }

	return _neighborGrid;
	// return _neighborGridTemp;
}