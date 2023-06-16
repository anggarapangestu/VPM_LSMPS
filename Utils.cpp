#include "Utils.hpp"

// Start the simulation counter number
void simUtil::startCounter(int step){
	this->iterDigitLen = 1;
	if (step != 0){
		iterDigitLen = 1 + std::floor(std::log10(step));
	}
	return;
}

// Printing the iteration number header
void simUtil::printHeader(int step){
	// Printing the iteration HEADER
	{
		printf("\n\n+");
		for(int _i = 0; _i < 16 - iterDigitLen/2; _i ++){
			printf("-");
		}
		printf(" Iteration Step %d ", step);
		for(int _i = 0; _i < 16 - iterDigitLen/2 - iterDigitLen%2; _i ++){
			printf("-");
		}
		printf("+\n");
	}
	return;
}

// Write the iteration saving file number
void simUtil::saveName(std::string & DataName, int step){
	int addDigit = Pars::maxDigitLen - this->iterDigitLen;
	for (int _spc = 0; _spc < addDigit; _spc++)
		DataName.append("0");
	DataName.append(std::to_string(step));
	return;
}

// Display the stability
void simUtil::stabilityEval(const Particle&par, std::vector<double>&max){
	// Initialize parameter
	double _courNum [3] = {0,0,0};  // [Current, Accumulative, Maximum]
	double _diffNum [3] = {0,0,0};  // [Current, Accumulative, Maximum]
	double _stabCrt [3] = {0,0,0};  // [Current, Accumulative, Maximum]
	for (int _i = 0; _i < par.num; _i++){
		// Calculating courant number
		_courNum[0] = std::sqrt(std::pow(par.u[_i],2) + std::pow(par.v[_i],2)) * Pars::dt / par.s[_i];
		_courNum[1] += _courNum[0];
		_courNum[2] = _courNum[2] > _courNum[0] ? _courNum[2] : _courNum[0];

		// Calculating diffusion number
		_diffNum[0] = Pars::NU * Pars::dt / (par.s[_i] * par.s[_i]);
		_diffNum[1] += _diffNum[0];
		_diffNum[2] = _diffNum[2] > _diffNum[0] ? _diffNum[2] : _diffNum[0];
		
		// Calculating stability criteria
		_stabCrt[0] = (std::abs(par.gz[_i]))/Pars::NU;
		_stabCrt[1] += _stabCrt[0];
		_stabCrt[2] = _stabCrt[2] > _stabCrt[0] ? _stabCrt[2] : _stabCrt[0];
	}
	// Average value
	_courNum[0] = _courNum[1] / par.num;
	_diffNum[0] = _diffNum[1] / par.num;
	_stabCrt[0] = _stabCrt[1] / par.num;
	
	// Displaying the value
	printf("Average courant number (C_av)           : %8.4f \n", _courNum[0]);
	printf("Average diffusion number (Phi_av)       : %8.4f \n", _diffNum[0]);
	printf("Average stability criteria (Re_h)       : %8.4f \n", _stabCrt[0]);
	printf("Max courant number (C_max)              : %8.4f \n", _courNum[2]);
	printf("Max diffusion number (Phi_max)          : %8.4f \n", _diffNum[2]);
	printf("Max stability criteria (Re_h)           : %8.4f \n", _stabCrt[2]);

	// Input the maximum value
	max.push_back(_courNum[2]);
	max.push_back(_diffNum[2]);
	max.push_back(_stabCrt[2]);
	return;
}

// Display the computational prdiction time
void simUtil::predictCompTime(int step, double curr_comp_time){
	// The value of "curr_comp_time" is in second (s)
	// Internal variable
	int est_time_d, est_time_h, est_time_m; double est_time_s;
	est_time_s = curr_comp_time * double(Pars::nt - step - 1);
	
	// Calculate Day
	est_time_d = int(est_time_s / (24 * 60 * 60));
	est_time_s -= est_time_d * (24 * 60 * 60);
	// Calculate Hour
	est_time_h = int(est_time_s / (60 * 60));
	est_time_s -= est_time_h * (60 * 60);
	// Calculate Minute
	est_time_m = int(est_time_s / (60));
	est_time_s -= est_time_m * (60);
	
	printf("\n<!> Estimation time to finish run:     %9.3f s", curr_comp_time*double(Pars::nt - step));
	if (est_time_d == 0){
		printf("\n<!> Estimation time to finish run: %2dh %2dm %5.2f s", est_time_h, est_time_m, est_time_s);
	}else{
		printf("\n<!> Estimation time to finish: %2dd %2dh %2dm %5.2f s", est_time_d, est_time_h, est_time_m, est_time_s);
	}
	return;
}


// The class to store the particle inside the cell
/*
std::vector<int> CellList::neighborCell(int pos){
	std::vector<int> _neighborList;
	
	// 2 dimension
	if (Pars::DIM == 2){
		int _botPos = pos - this->nx;
		int _topPos = pos + this->nx;
		int _operatorPos[8] = {_botPos - 1, _botPos, _botPos + 1,
							pos - 1, pos+ 1,
							_topPos - 1, _topPos, _topPos + 1};
		for (int i = 0; i < 8; i++){
			if (_operatorPos[i] < 0)
			_neighborList.emplace_back(0);
		}
	}
	// 3 dimension
	else if (Pars::DIM == 3){
		for (int i = 0; i < 26; i++){
			_neighborList.emplace_back(0);
		}
	}
}

// The class to store the particle inside the cell
void CellList::createCellList(const Particle &par){
	// Internal variable
	double *pos_min = new double[3];
	double *pos_max = new double[3];
	double xmin, ymin, zmin;

	// Label
	this->h = Pars::sigma * Pars::r_sup + Pars::U_inf * Pars::dt * 1.2;

	// find boundary position of each element in particle
	pos_min[0] = par.x[0];
	pos_min[1] = par.y[0];
	// pos_min[2] = par.z[0];
	pos_max[0] = par.x[0];
	pos_max[1] = par.y[0];
	// pos_max[2] = par.z[0];
	for (int i = 1; i < par.num; i++){
		pos_min[0] = pos_min[0] < par.x[i] ? pos_min[0] : par.x[i];
		pos_min[1] = pos_min[1] < par.y[i] ? pos_min[1] : par.y[i];
		pos_min[2] = pos_min[2] < par.z[i] ? pos_min[2] : par.z[i];
		pos_max[0] = pos_max[0] > par.x[i] ? pos_max[0] : par.x[i];
		pos_max[1] = pos_max[1] > par.y[i] ? pos_max[1] : par.y[i];
		pos_max[2] = pos_max[2] > par.z[i] ? pos_max[2] : par.z[i];
	}
	
	this->nx = std::ceil((20 * Pars::sigma + pos_max[0] - pos_min[0])/this->h);     // Number of cell in x direction
	this->ny = std::ceil((pos_max[1] - pos_min[1])/this->h);     // Number of cell in y direction
	this->nz = std::ceil((pos_max[2] - pos_min[2])/this->h);     // Number of cell in z direction

	// xmin = - Pars::xdom - 10 * Pars::sigma
	// xmax = - Pars::xdom - 10 * Pars::sigma + this->nx*this->h
	// ymin = - this->ny*this->h / 2.0
	// ymax = this->ny*this->h / 2.0
	// zmin = - this->nz*this->h / 2.0
	// zmax = this->nz*this->h / 2.0
	// this->num = this->nx * this->ny;    // * this->nz
	
	// Put the particle inside the cell
	int pos;
	std::vector<std::vector<int>> child;;
	this->particle_inside.resize(this->num, std::vector<int>(0));
	for (int i = 0; i < par.num; i++){
		// Check the position of the particle then put inside the box
		double x = par.x[i] + Pars::xdom + 10 * Pars::sigma;
		double y = par.y[i] + this->ny*this->h / 2.0;
		double z = par.z[i] + this->nz*this->h / 2.0;
		pos = std::floor(x / this->h) + this->nx * std::floor(y / this->h) + this->nx * this->ny * std::floor(z / this->h);
		std::vector<std::vector<int>> child;
    	std::vector<std::vector<int>> particle_inside;
		
		particle_inside[pos].emplace_back(i);
	}
}
*/