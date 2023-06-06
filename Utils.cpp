#include "Utils.hpp"
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

// // constructor
// // Body::Body();
// // methods
// // 	 getter
// const std::vector<double>& Body::get_x() {return _x; }
// const std::vector<double>& Body::get_y() {return _y; }
// const std::vector<double>& Body::get_uT() {return _uT; }
// const std::vector<double>& Body::get_vT() {return _vT; }
// const std::vector<double>& Body::get_uR() {return _uR; }
// const std::vector<double>& Body::get_vR() {return _vR; }
// const std::vector<double>& Body::get_uDEF() {return _uDEF; }
// const std::vector<double>& Body::get_vDEF() {return _vDEF; }
// // 	setter
// void Body::set_x(const std::vector<double> &x) {_x = x; }
// void Body::set_y(int size) {_y.clear(); _y.resize(size); }
// void Body::set_uT(int size) {_uT.clear(); _uT.resize(size); }
// void Body::set_vT(int size) {_vT.clear(); _vT.resize(size); }
// void Body::set_uR(int size) {_uR.clear(); _uR.resize(size); }
// void Body::set_vR(int size) {_vR.clear(); _vR.resize(size); }
// void Body::set_uDEF(int size) {_uDEF.clear(); _uDEF.resize(size); }
// void Body::set_vDEF(int size) {_vDEF.clear(); _vDEF.resize(size); }
// // destructor
// // Body::~Body();
*/