#include "pressure_poisson.hpp"

#ifndef INCLUDED_REMESHING_BASE
#include "../remeshing/base_remeshing.hpp"
#endif

void pressure_poisson::get_pressure(Particle& p){
	// Update the start clock
	clock_t _time = clock();
				
	// Solver computation
	printf("\nCalculate Pressure ...\n");
	
	// Resize the velocity of particle
	p.P.clear();
	p.P.resize(p.num,0.0);

	// ******************************* NEW FMM DATA *******************************
	// Internal variable
	fmm2D FMM_step;			// The FMM operator class
	std::vector<std::vector<double>> particle_POS(p.num, std::vector<double>(Pars::DIM,0.0));	// Particle position
	std::vector<double> particle_SRC(p.num, 0.0);	// The particle source
	std::vector<bool> particle_mark(p.num, false);	// The particle active mark
	double start,finish;							// Computational time
	
	start = omp_get_wtime();
	// [1] Calculcate the pressure source
    std::vector<double> _SRC(p.num, 0.0);
	LSMPSa lsmpsa_du;     // module to calculate spatial differential of u velocity
	LSMPSa lsmpsa_dv;     // module to calculate spatial differential of v velocity

	lsmpsa_du.set_LSMPS(p.x, p.y, p.s, p.u, p.neighbor);
	lsmpsa_dv.set_LSMPS(p.x, p.y, p.s, p.v, p.neighbor);
	std::vector<double> _dudx = lsmpsa_du.get_ddx();
	std::vector<double> _dudy = lsmpsa_du.get_ddy();
	std::vector<double> _dvdx = lsmpsa_dv.get_ddx();

    // Calculcate the divergenve pressure source
    for (size_t i = 0; i < p.num; i++){
        // Here calculate 2 * (dudx^2 + dudy * dvdx)
        _SRC[i] = 2.0 * (_dudx[i] * _dudx[i] + _dudy[i] * _dvdx[i]);
    }

	// Determine active source
	// Here find the maximum value
	double SRC_max;
	SRC_max = std::abs(_SRC[0]);
    for (size_t i = 1; i < p.num; i++){
        SRC_max = SRC_max > std::abs(_SRC[i]) ? SRC_max : std::abs(_SRC[i]);
    }

	// ===================================================
    // Assign each variable
	for (size_t i = 0; i < p.num; i++){
		// Particle position
		particle_POS[i][0] = p.x[i];
		particle_POS[i][1] = p.y[i];

		// Particle source
		particle_SRC[i] = -(Pars::RHO / (2 * M_PI)) * _SRC[i] * (p.s[i] * p.s[i]);
		
		// Particle mark
		particle_mark[i] = std::abs(_SRC[i]) > SRC_max * std::pow(10,-4) ? true : false;
		particle_mark[i] = particle_mark[i] | p.isActive[i];
	}

	// // DEBUGGING: save data value
	// this->save_pressure_source(p.x,p.y,_SRC,particle_mark);

	finish = omp_get_wtime();
	printf("<+> Initialization   : %f s\n", finish-start);

	// Initialize the tree
	start = omp_get_wtime();
	
	// Initialization Tree Map
	// treeData.initializeTree(particle_POS);
	// treeData.createTree(particle_POS, particle_mark);
	// if (step == 0){
		treeData.initializeTree(particle_POS);
		treeData.createTree(particle_POS, particle_mark);
	// }else{
	// 	treeData.updateCell(particle_POS, particle_mark);
	// }

	// Tree Map Saved
	if (false){
		std::string name = "P";
		name += std::to_string(this->counter);
		this->counter ++;
		treeData.saveTreeCell(name);	
	}

	finish = omp_get_wtime();
	printf("<+> Tree finished in : %f s\n", finish-start);

	// Calculate the tree data
	FMM_step.calcPotential(treeData, particle_POS, particle_mark, particle_SRC);
	
	// Get the final result
	std::vector<double> P;
	FMM_step.get_Potential(P);

	// ===================================================
	// Update the velocity
	for (size_t i = 0; i < p.num; i++)
	{
		p.P[i] = P[i];
		// p.P[i] = _SRC[i];
		
	}

	// ******************************* NEW FMM DATA *******************************

	// Display computational time
	_time = clock() - _time;    // Calculate the time duration
	printf("<-> Velocity calculation comp. time:   [%f s]\n", (double)_time/CLOCKS_PER_SEC);

	
	// ******************************* Boundary Limit *******************************
	// ******************************* Boundary Limit *******************************
	// Calculation of the pressure coefficient
	int _num = 76;
	std::vector<double> alpha(_num,0.0);
	Particle surface, buffer;
	surface.num = _num;
	surface.x.resize(_num,0.0);
	surface.y.resize(_num,0.0);
	surface.s.resize(_num,Pars::sigma);
	surface.P.resize(_num,0.0);
	for (int i = 0; i < _num; i++){
		surface.x[i] = -(0.5+0*Pars::sigma)*std::cos((i*2.0*M_PI)/_num);
		surface.y[i] = (0.5+0*Pars::sigma)*std::sin((i*2.0*M_PI)/_num);
	}
	// Definition of new set particle
	// Calculation of the pressure coefficient
	buffer.num = 0;
	buffer.x.clear();
	buffer.y.clear();
	buffer.s.clear();
	buffer.P.clear();
	double lim_max = 0.5 + Pars::sigma * 7;
	double lim_min = -0.5 - Pars::sigma * 7;
	for (int i = 0; i < p.num; i++){
		if ((p.x[i] > lim_min && p.x[i] < lim_max) && (p.y[i] > lim_min && p.y[i] < lim_max)){
			buffer.num++;
			buffer.x.push_back(p.x[i]);
			buffer.y.push_back(p.y[i]);
			buffer.s.push_back(p.s[i]);
			buffer.P.push_back(p.P[i]);
		}
	}
	
	// Find the neighbor for calculation use direct search
	std::vector<std::vector<int>> _par2grdNeighbor,_grd2grdNeighbor;
	base_remeshing d_base_remeshing;
	_par2grdNeighbor = d_base_remeshing.inter_search(buffer, surface);
	_grd2grdNeighbor = d_base_remeshing.inter_search(surface, surface);

	// Interpolation of vector
	LSMPSb lsmps_p;
	lsmps_p.set_LSMPS(surface.x, surface.y, surface.s, surface.P,
                     buffer.x, buffer.y, buffer.s, buffer.P,
                     _par2grdNeighbor, _grd2grdNeighbor);
	surface.P = lsmps_p.get_d00();

	// Save the data
	std::ofstream _data;
	std::string name = "output/Cp_";
	name += std::to_string(this->counter1);
	name += ".csv";
	this->counter1 ++;
	_data.open(name);
	_data << "x,y,Cp\n";
	
	for (size_t i = 0; i < surface.num; i++){
		_data << surface.x[i] << ","
		      << surface.y[i] << ","
			  << surface.P[i] << "\n";
	}
	_data.close();
	return;
}

void pressure_poisson::save_pressure_source(const std::vector<double>&x, const std::vector<double>&y, const std::vector<double>&src, const std::vector<bool>&mark){
	std::ofstream _data;
	std::string name = "output/pressureSRC";
	name += std::to_string(this->counter);
	name += ".csv";
	counter ++;
	_data.open(name);
	_data << "x,y,SRC,mark\n";
	
	for (size_t i = 0; i < x.size(); i++){
		_data << x[i] << ","
		      << y[i] << ","
			  << src[i] << ","
			  << mark[i] << "\n";
	}
	_data.close();
	return;
}