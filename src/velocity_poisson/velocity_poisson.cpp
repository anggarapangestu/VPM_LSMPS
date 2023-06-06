#include "velocity_poisson.hpp"

// =================================================================
//  TODO: solving POISSON EQUATION for Particles by using Biot_Savart and technique FMM to accelerate computation speed
//  (VIC : solves Poisson Equation for grid-nodes by using FDM-finite differental method and technique FFM- Fast Furier Transform to accelerate computation speed)
//  input: PARTICLE info: vortecity-strength, position
//  output: calculate velocity and Stretching term dgdt_st_x
// =================================================================

void velocity_poisson::get_velocity(Particle &p, const int step)
{
	// Update the start clock
	clock_t _time = clock();
				
	// Solver computation
	printf("\nCalculate Velocity  ...\n");
	printf("<+> Calculate the rotational velocity term\n");
	
	// Resize the velocity of particle
	p.u.clear();
	p.u.resize(p.num,0.0e0);
	p.v.clear();
	p.v.resize(p.num,0.0e0);
	
	// 	// -- internal variables
	// 	const static int n_inter = 1; // ! change into class static variables
	// 	std::vector<int> _index;
	// 	Particle _particle;		  // store current active particles box (to increase robustness)
	// 	Particle _particleDense;  // store particle with core size as smallest from current active particles
	// 	_particle.num = 0;
	// 	_particleDense.num = 0;
	// 	velocity_biot_savart FMM;

		
	// 	// TODO: store current active particles (Ngga efficient)
	// 	for (size_t i = 0; i < p.num; i++)
	// 	{
	// 		//if (p.isActive[i] == true)
	// 		//{
	// 			_particle.x.push_back(p.x[i]);
	// 			_particle.y.push_back(p.y[i]);
	// 			_particle.s.push_back(p.s[i]);
	// 			_particle.gz.push_back(p.gz[i]);// * std::pow(p.s[i] / Pars::sigma, 2));
	// 			_particle.u.push_back(0.0e0);
	// 			_particle.v.push_back(0.0e0);
	// 			_index.push_back(i);
	// 			_particle.num++;
	// 			_particleDense.x.push_back(p.x[i]);
	// 			_particleDense.y.push_back(p.y[i]);
	// 			_particleDense.s.push_back(p.s[i]);
	// 			_particleDense.gz.push_back(p.gz[i]);// * std::pow(p.s[i] / Pars::sigma, 2));
	// 			_particleDense.u.push_back(p.u[i]);
	// 			_particleDense.v.push_back(p.v[i]);
	// 			_particleDense.num++;
	// 		//}

	// 	}


	// 	// TODO: calculate FMM by calling Fortran subroutine
	// 	// := begin
	// 	int np = _particle.num;
	// 	double *x = new double[_particle.num];
	// 	double *y = new double[_particle.num];
	// 	double *s = new double[_particle.num];
	// 	double *gz = new double[_particle.num];
	// 	double *u = new double[_particle.num];
	// 	double *v = new double[_particle.num];
	// 	/*
	// 	for (int i = 0; i < _particle.num; i++)
	// 	{
	// 		x[i] = _particle.x[i];
	// 		y[i] = _particle.y[i];
	// 		s[i] = _particle.s[i];
	// 		gz[i] = _particle.gz[i];
	// 		u[i] = 0.0e0;
	// 		v[i] = 0.0e0;
	// 	}*/

	// 	int _n0 = 1;
	// 	int _iCutoff = Pars::icutoff;
	// 	int _nS = Pars::n_s;
	// 	int _nInter = n_inter;
	// 	int _ndp = Pars::ndp;

	// 	// printf("<+> Called from Fortran90 Utility\n");
	// 	/*FortranUtils::biotsavart_fmm_2d_(&_n0, &np, &np, x, y, s, u, v,
	// 									&_n0, &np, &np, x, y, s, gz,
	// 									&_iCutoff, &_nS, &_nInter, &_ndp);
	// 	for (int i = 0; i < _particle.num; p)
	// 	{
	// 		_particle.u[i] = u[i];
	// 		_particle.v[i] = v[i];
	// 	}*/
		
	// 	FMM.biotsavart_fmm_2d( _particle, _particleDense,  _iCutoff, _nS,  _nInter,  _ndp);
	// 	//FMM.biotsavart_direct_2d(_particle, _particleDense);

	// 	// update base particle velocity
	// 	// Dapetnya kecepatan
	// 	for (size_t i = 0; i < _particle.num; i++)
	// 	{
	// 		int _isActiveIndex = _index[i];
	// 		p.u[_isActiveIndex] = _particle.u[i];
	// 		p.v[_isActiveIndex] = _particle.v[i];
	// 	}
	// delete[] x;
	// delete[] y;
	// delete[] s;
	// delete[] gz;
	// delete[] u;
	// delete[] v;

	// ******************************* NEW FMM DATA *******************************
	// Internal variable
	fmm2D FMM_step;			// The FMM operator class
	std::vector<std::vector<double>> particle_POS(p.num, std::vector<double>(Pars::DIM,0.0));	// Particle position
	std::vector<double> particle_SRC(p.num, 0.0);	// The particle source
	std::vector<bool> particle_mark(p.num, false);	// The particle active mark
	double start,finish,total;						// Computational time
	
	// Assign each variable
	for (size_t i = 0; i < p.num; i++){
		// Particle position
		particle_POS[i][0] = p.x[i];
		particle_POS[i][1] = p.y[i];

		// Particle source
		particle_SRC[i] = - p.gz[i] / (2 * M_PI);
		
		// Particle mark
		particle_mark[i] = p.isActive[i];
	}

	// Initialize the tree
	start = omp_get_wtime();
	
	// Initialization Tree Map
	// treeData.initializeTree(particle_POS);
	// treeData.createTree(particle_POS, particle_mark);
	if (step == 0){
		treeData.initializeTree(particle_POS);
		treeData.createTree(particle_POS, particle_mark);
	}else{
		treeData.updateCell(particle_POS, particle_mark);
	}

	// Tree Map Saved
	if (false){
		std::string name = "V";
		name += std::to_string(step);
		treeData.saveTreeCell(name);	
	}

	finish = omp_get_wtime();
	printf("<+> Tree finished in : %f s\n", finish-start);

	// Calculate the tree data
	FMM_step.calcField(treeData, particle_POS, particle_mark, particle_SRC);

	// // Calculate the tree data
	// FMM_step.calcPotential(treeData, particle_POS, particle_mark, particle_SRC);
	
	// Get the final result
	std::vector<double> Ex;
	std::vector<double> Ey;
	FMM_step.get_Field(Ex, Ey);

	// Update the velocity
	for (size_t i = 0; i < p.num; i++)
	{
		p.u[i] = Ey[i];
		p.v[i] = -Ex[i];
	}

	// ******************************* NEW FMM DATA *******************************

	// Display computational time
	_time = clock() - _time;    // Calculate the time duration
	printf("<-> Velocity calculation comp. time:   [%f s]\n", (double)_time/CLOCKS_PER_SEC);

	// := end

	// double *_maxPar = new double[2];
	// double *_minPar = new double[2];
	// _minPar[0] = *std::min_element(_particleActive.x.begin(), _particleActive.x.end()) - 10 * Pars::sigma;
	// _minPar[1] = *std::min_element(_particleActive.y.begin(), _particleActive.y.end()) - 10 * Pars::sigma;
	// _maxPar[0] = *std::max_element(_particleActive.x.begin(), _particleActive.x.end()) + 10 * Pars::sigma;
	// _maxPar[1] = *std::max_element(_particleActive.y.begin(), _particleActive.y.end()) + 10 * Pars::sigma;
	// for (size_t i = 0; i < p.num; i++)
	// {
	// 	if (p.x[i] >= _minPar[0] && p.x[i] <= _maxPar[0] &&
	// 		p.y[i] >= _minPar[1] && p.y[i] <= _maxPar[1])
	// 	{
	// 		_particle.x.push_back(p.x[i]);
	// 		_particle.y.push_back(p.y[i]);
	// 		_particle.s.push_back(p.s[i]);
	// 		_particle.gz.push_back(p.gz[i] * std::pow(p.s[i] / Pars::sigma, 2));
	// 		_particle.u.push_back(p.u[i]);
	// 		_particle.v.push_back(p.v[i]);
	// 		_index.push_back(i);
	// 		_particle.num++;
	// 	}
	// }
	// delete[] _maxPar;
	// delete[] _minPar;

	// TODO: create base mesh as lowest level to be calculated using FMM
	// procedures:
	// Step 1: generate lowest level mesh [lv1]
	// Step 2: interpolate active particles [lv0] properties (i.e. circulation) into lv1
	// Step 3: calculate Poisson equation from lv1
	// Step 4: interpolate velocity from lv1 back to lv0

	// !STEP 1: begin
	// double *_maxDom = new double[2];
	// double *_minDom = new double[2];
	// int *_nNode = new int[2];
	// _minDom[0] = *std::min_element(_particle.x.begin(), _particle.x.end()) - 10 * Pars::sigma;
	// _minDom[1] = *std::min_element(_particle.y.begin(), _particle.y.end()) - 10 * Pars::sigma;
	// _maxDom[0] = *std::max_element(_particle.x.begin(), _particle.x.end()) + 10 * Pars::sigma;
	// _maxDom[1] = *std::max_element(_particle.y.begin(), _particle.y.end()) + 10 * Pars::sigma;

	// double _ratio = 1.0e0;
	// _nNode[0] = std::ceil((_maxDom[0] - _minDom[0]) / (_ratio * Pars::sigma));
	// _nNode[1] = std::ceil((_maxDom[1] - _minDom[1]) / (_ratio * Pars::sigma));
	// // assign grid
	// for (size_t i = 0; i < _nNode[0]; i++)
	// {
	// 	double _x = _minDom[0] + static_cast<double>(i) * _ratio * Pars::sigma;
	// 	for (size_t j = 0; j < _nNode[1]; j++)
	// 	{
	// 		double _y = _minDom[1] + static_cast<double>(j) * _ratio * Pars::sigma;
	// 		_particleDense.x.push_back(_x);
	// 		_particleDense.y.push_back(_y);
	// 		_particleDense.s.push_back(_ratio * Pars::sigma);
	// 		_particleDense.num++;
	// 	}
	// }
	// _particleDense.u.resize(_particleDense.num);
	// _particleDense.v.resize(_particleDense.num);
	// _particleDense.gz.resize(_particleDense.num);

	// delete[] _maxDom;
	// delete[] _minDom;
	// delete[] _nNode;
	// // generate nighborhood
	// // current notation is:
	// // 		active particle <- particle
	// // 		lowest level mesh <- grid
	// // neighborhood consist of:
	// // 1. particle--particle
	// // 2. grid--particle
	// // 3. particle--grid
	// // 4. grid--grid
	// std::vector<std::vector<int>> _par2parNeighbor = _base_remeshing.inter_search(_particle, _particle);		   // 2
	// std::vector<std::vector<int>> _grd2parNeighbor = _base_remeshing.inter_search(_particleDense, _particle);	  // 2
	// std::vector<std::vector<int>> _par2grdNeighbor = _base_remeshing.inter_search(_particle, _particleDense);	  // 1
	// std::vector<std::vector<int>> _grd2grdNeighbor = _base_remeshing.inter_search(_particleDense, _particleDense); // 1
	// !STEP 1: end

	// !STEP 2: begin
	// _lsmpsb.set_LSMPS(
	// 	_particleDense.x, _particleDense.y, _particleDense.s, _particleDense.gz,
	// 	_particle.x, _particle.y, _particle.s, _particle.gz,
	// 	_par2grdNeighbor, _grd2grdNeighbor);
	// _particleDense.gz = _lsmpsb.get_d00();
	// !STEP 2: end

	// printf("Solving Poisson's equation ... \n");
	// // !STEP 3: begin
	// if (Pars::ioptfmm == 0)
	// {
	// 	// biot_savart.biotsavart_direct_2d(p, p);
	// 	biot_savart.biotsavart_direct_2d(_particle, _particle);
	// 	// biot_savart.biotsavart_direct_2d(_particleDense, _particleDense);
	// }
	// else if (Pars::ioptfmm == 1)
	// {
	// 	// biot_savart.biotsavart_fmm_2d(p, p, Pars::icutoff, Pars::n_s, n_inter, Pars::ndp);
	// 	biot_savart.biotsavart_fmm_2d(_particle, _particle, Pars::icutoff, Pars::n_s, n_inter, Pars::ndp);
	// 	// biot_savart.biotsavart_fmm_2d(_particleDense, _particleDense, Pars::icutoff, Pars::n_s, n_inter, Pars::ndp);
	// }
	// !STEP 3: end

	// !STEP 4: begin
	// _lsmpsb.set_LSMPS(
	// 	_particle.x, _particle.y, _particle.s, _particle.u,
	// 	_particleDense.x, _particleDense.y, _particleDense.s, _particleDense.u,
	// 	_grd2parNeighbor, _par2parNeighbor);
	// _particle.u = _lsmpsb.get_d00();

	// _lsmpsb.set_LSMPS(
	// 	_particle.x, _particle.y, _particle.s, _particle.v,
	// 	_particleDense.x, _particleDense.y, _particleDense.s, _particleDense.v,
	// 	_grd2parNeighbor, _par2parNeighbor);
	// _particle.v = _lsmpsb.get_d00();

	// // update base particle velocity
	// // Dapetnya kecepatan
	// for (size_t i = 0; i < _particle.num; i++)
	// {
	// 	int _isActiveIndex = _index[i];
	// 	p.u[_isActiveIndex] = _particle.u[i];
	// 	p.v[_isActiveIndex] = _particle.v[i];
	// }
	// !STEP 4: end
}
