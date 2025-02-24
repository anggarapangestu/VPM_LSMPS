#include "stability.hpp"
#include "../LSMPS/LSMPSa.hpp"

/**
 *  @brief  Display the stability into the console.
 *         
 *  @param  _particle The particle data for stability calculaltion.
 *  @param  _maxValue [OUTPUT] Stability criteria global maximum value.
 *  @param  _step  The current iteration step
*/
void stability::stabilityEval(const Particle &par, std::vector<double*> &max, const int &step){
	// A prompt in stability evaluation

    // Procedure
    // [1] Calculate the stability value at current iteration
    // [2] Update the moving time integration value
    // [3] Save stability data

	printf("Evaluating stability ...\n");
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

	// Internal container
	// [(1)Current Evaluated -> Average, (2)Accumulative in Domain, (3)Global Maximum]
	double _courNum [3]    = {0,0,0};  // The courant number container
	double _diffNum [3]    = {0,0,0};  // The diffusion number container
	double _vorMeshNum [3] = {0,0,0};  // The vortex mesh reynold number container
    double _lagStrNum [3]  = {0,0,0};  // The lagrangian stretching criteria container

    // Calculate the evaluation
    if(Pars::stab_courant) this->courant_eval(par, _courNum);
    if(Pars::stab_diffusion) this->diffusion_eval(par, _diffNum);
    // if(Pars::stab_vortex) this->vortex_eval(par, _vorMeshNum, 1);		// Type 1: Evaluate the combined diffusion and lagrangian
    if(Pars::stab_vortex) this->vortex_eval(par, _vorMeshNum, 2);			// Type 2: Evaluate the langrangian by vorticity
    if(Pars::stab_lag_stretch) this->lag_str_eval(par, _lagStrNum);

	// Determine the current step is initialization or not
	bool init = false;
	if      (Pars::opt_start_state == 0) init = (step == 0);
	else if (Pars::opt_start_state == 1) init = (step == (Pars::resume_step + 1));

	// Update the global stability criteria maximum value
	if (init == true){
		*(max[0]) = Pars::stab_courant ? _courNum[2] : 0.0;
		*(max[1]) = Pars::stab_diffusion ? _diffNum[2] : 0.0;
		*(max[2]) = Pars::stab_vortex ? _vorMeshNum[2] : 0.0;
        *(max[3]) = Pars::stab_lag_stretch ? _lagStrNum[2] : 0.0;
	}else{
		*(max[0]) = Pars::stab_courant ? std::max(*(max[0]), _courNum[2]) : 0.0;
		*(max[1]) = Pars::stab_diffusion ? std::max(*(max[1]), _diffNum[2]) : 0.0;
		*(max[2]) = Pars::stab_vortex ? std::max(*(max[2]), _vorMeshNum[2]) : 0.0;
        *(max[3]) = Pars::stab_lag_stretch ? std::max(*(max[3]), _lagStrNum[2]) : 0.0;
	}

	// Save the stability data
	if (Pars::flag_save_stability == true){
		// Initialize the writter
		std::ofstream _write;
		std::string fileName = "output/Stability_History";
		if      (Pars::opt_start_state == 0) fileName += ".csv";
		else if (Pars::opt_start_state == 1) fileName += "_Resume.csv";
		
		// Save the header
		if (init == true){
			_write.open(fileName);
			_write << "time,C,Phi,Re_h,C_L\n";
			_write.close();
		}
		
		// Save the data
		_write.open(fileName, std::ofstream::out | std::ofstream::app);
		_write <<  "" << step*Pars::dt 
			   << "," << _courNum[2]
			   << "," << _diffNum[2]
			   << "," << _vorMeshNum[2]
               << "," << _lagStrNum[2]
			   << "\n";
		_write.close();
	}

	// Particle generation initialization summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Stability evaluation comp. time.   [%f s]\n\n", _time);

	// Displaying the value
	if(Pars::stab_courant)      printf("Average courant number (C_av)           : %8.4f \n", _courNum[0]);
	if(Pars::stab_diffusion)    printf("Average diffusion number (Phi_av)       : %8.4f \n", _diffNum[0]);
	if(Pars::stab_vortex)       printf("Average vortex mesh criteria (Re_h)     : %8.4f \n", _vorMeshNum[0]);
    if(Pars::stab_lag_stretch)  printf("Average lag. stretch criteria (C_L)     : %8.4f \n", _lagStrNum[0]);
    
    printf("---------------------------------------------------\n");
	if(Pars::stab_courant)      printf("Max courant number (C_max) [< 1.0]      : %8.4f \n", _courNum[2]);
	if(Pars::stab_diffusion)    printf("Max diffusion number (Phi_max) [< 0.5]  : %8.4f \n", _diffNum[2]);
	if(Pars::stab_vortex)       printf("Max vortex mesh criteria (Re_h) [O(1)]  : %8.4f \n", _vorMeshNum[2]);
    if(Pars::stab_lag_stretch)  printf("Max lag. stretch criteria (C_L) [< 0.5] : %8.4f \n", _lagStrNum[2]);

	printf("---------------------------------------------------\n");
	if(Pars::stab_courant)      printf("Lifetime Courant (C_max)    	[< 1.0] : %8.4f \n", *(max[0]));
	if(Pars::stab_diffusion)    printf("Lifetime Diffusion (Phi_max)  	[< 0.5] : %8.4f \n", *(max[1]));
	if(Pars::stab_vortex)       printf("Lifetime Vortex Mesh (Re_h) 	[O(1)]  : %8.4f \n", *(max[2]));
    if(Pars::stab_lag_stretch)  printf("Lifetime Lag. Stretch (C_L_max) [< 0.5] : %8.4f \n", *(max[3]));
	
	return;
}


/**
 *  @brief The evaluation calculation of courant number (C)
 *  @param par The particle to be evaluated
 *  @param _courNum The courant number container
 */
void stability::courant_eval(const Particle par, double* _courNum){    
    // Calculation operation
	#if (DIM == 2)
		for (int _i = 0; _i < par.num; _i++){
			// Calculating courant number
			double Velocity = std::sqrt(par.u[_i]*par.u[_i] + par.v[_i]*par.v[_i]);
			_courNum[0] = Velocity * Pars::dt / par.s[_i];
			_courNum[1] += _courNum[0];
			_courNum[2] = _courNum[2] > _courNum[0] ? _courNum[2] : _courNum[0];
		}
	#elif (DIM == 3)
		for (int _i = 0; _i < par.num; _i++){
			// Calculating courant number
			double Velocity = std::sqrt(par.u[_i]*par.u[_i] + par.v[_i]*par.v[_i] + par.w[_i]*par.w[_i]);
			_courNum[0] = Velocity * Pars::dt / par.s[_i];
			_courNum[1] += _courNum[0];
			_courNum[2] = _courNum[2] > _courNum[0] ? _courNum[2] : _courNum[0];
		}
	#endif

    // Calculate the average value (put into the current evaluation container)
	_courNum[0] = _courNum[1] / par.num;
    return;
}

/**
 *  @brief The evaluation calculation of diffusion number (Phi)
 *  @param par The particle to be evaluated
 *  @param _diffNum The diffusion number container
 */
void stability::diffusion_eval(const Particle par, double* _diffNum){    
    // Calculation operation
	#if (DIM == 2)
		for (int _i = 0; _i < par.num; _i++){
			// Calculating diffusion number
			_diffNum[0] = Pars::NU * Pars::dt / (par.s[_i] * par.s[_i]);
			_diffNum[1] += _diffNum[0];
			_diffNum[2] = _diffNum[2] > _diffNum[0] ? _diffNum[2] : _diffNum[0];
		}
	#elif (DIM == 3)
		for (int _i = 0; _i < par.num; _i++){
			// Calculating diffusion number
			_diffNum[0] = Pars::NU * Pars::dt / (par.s[_i] * par.s[_i]);
			_diffNum[1] += _diffNum[0];
			_diffNum[2] = _diffNum[2] > _diffNum[0] ? _diffNum[2] : _diffNum[0];
		}
	#endif

    // Calculate the average value (put into the current evaluation container)
	_diffNum[0] = _diffNum[1] / par.num;
    return;
}

/**
 *  @brief The evaluation calculation of vortex number (ReH)
 *  @param par The particle to be evaluated
 *  @param _vrtMeshNum The vortex number container
 *  @param _type The type of vortex number calculation [1:= Calculate vortex strength, 2:= Calculate the vorticity]
 */
void stability::vortex_eval(const Particle par, double* _vrtMeshNum, int _type){
    if (_type == 1){
		// Calculation operation by vortex strength
		#if (DIM == 2)
			for (int _i = 0; _i < par.num; _i++){
				// Calculating vorticity reynolds criteria
				_vrtMeshNum[0] = (std::abs(par.gz[_i])) / Pars::NU;
				_vrtMeshNum[1] += _vrtMeshNum[0];
				_vrtMeshNum[2] = _vrtMeshNum[2] > _vrtMeshNum[0] ? _vrtMeshNum[2] : _vrtMeshNum[0];
			}
		#elif (DIM == 3)
			for (int _i = 0; _i < par.num; _i++){
				// Calculating stability criteria
				// double Vorticity = std::sqrt(par.vortx[_i]*par.vortx[_i] + par.vorty[_i]*par.vorty[_i] + par.vortz[_i]*par.vortz[_i]);
				_vrtMeshNum[0] = std::abs(par.vorticity[_i]) * par.s[_i] * par.s[_i] / Pars::NU;
				_vrtMeshNum[1] += _vrtMeshNum[0];
				_vrtMeshNum[2] = _vrtMeshNum[2] > _vrtMeshNum[0] ? _vrtMeshNum[2] : _vrtMeshNum[0];
			}
		#endif
	}
	
	else if(_type == 2){
		// Calculation operation by vorticity
		#if (DIM == 2)
			for (int _i = 0; _i < par.num; _i++){
				// Calculating vorticity reynolds criteria
				_vrtMeshNum[0] = std::abs(par.vorticity[_i]) * Pars::dt;
				_vrtMeshNum[1] += _vrtMeshNum[0];
				_vrtMeshNum[2] = _vrtMeshNum[2] > _vrtMeshNum[0] ? _vrtMeshNum[2] : _vrtMeshNum[0];
			}
		#elif (DIM == 3)
			for (int _i = 0; _i < par.num; _i++){
				_vrtMeshNum[0] = std::abs(par.vorticity[_i]) * Pars::dt;
				_vrtMeshNum[1] += _vrtMeshNum[0];
				_vrtMeshNum[2] = _vrtMeshNum[2] > _vrtMeshNum[0] ? _vrtMeshNum[2] : _vrtMeshNum[0];
			}
		#endif
	}

	else{
		throw std::runtime_error("Error: The type is outside the range!");
	}

    // Calculate the average value (put into the current evaluation container)
	_vrtMeshNum[0] = _vrtMeshNum[1] / par.num;
    return;
}

/**
 *  @brief The evaluation calculation of lagrangian stretching criteria (C_L)
 *  @param _p The particle to be evaluated
 *  @param _lagStrNum The lagrangian stretching criteria container
 */
void stability::lag_str_eval(const Particle _p, double* _lagStrNum){	
	// Calculation operation
	#if (DIM == 2)
        // Calculate the velocity first order spatial differential
        std::vector<double> dudx, dudy, dvdx, dvdy;
        LSMPSa lsmpsVel2D;
        
        // Calculate the spatial differential for x-velocity 
        lsmpsVel2D.set_LSMPS(_p.x, _p.y, _p.s, _p.u, 
                           _p.x, _p.y, _p.s, _p.u, _p.neighbor);
        dudx = lsmpsVel2D.get_ddx();
        dudy = lsmpsVel2D.get_ddy();

        // Calculate the spatial differential for y-velocity 
        lsmpsVel2D.set_LSMPS(_p.x, _p.y, _p.s, _p.v, 
                           _p.x, _p.y, _p.s, _p.v, _p.neighbor);
        dvdx = lsmpsVel2D.get_ddx();
        dvdy = lsmpsVel2D.get_ddy();

		// Iteration for the velocity value
        for (int _i = 0; _i < _p.num; _i++){
			std::vector<double> velGrad(4);
            velGrad[0] = std::abs(dudx[_i]);
            velGrad[1] = std::abs(dudy[_i]);
            velGrad[2] = std::abs(dvdx[_i]);
            velGrad[3] = std::abs(dvdy[_i]);

            // Find the maximum gradient value
            double maxGrad = *std::max_element(velGrad.begin(), velGrad.end());
			
			_lagStrNum[0] = maxGrad * Pars::dt;
			_lagStrNum[1] += _lagStrNum[0];
			_lagStrNum[2] = _lagStrNum[2] > _lagStrNum[0] ? _lagStrNum[2] : _lagStrNum[0];
		}

	#elif (DIM == 3)
        std::vector<double> dudx, dudy, dudz;       // The differential on x-velocity
        std::vector<double> dvdx, dvdy, dvdz;       // The differential on y-velocity
        std::vector<double> dwdx, dwdy, dwdz;       // The differential on z-velocity
        LSMPSa lsmpsVel3D;
        
        // Calculate the spatial differential for x-velocity 
        lsmpsVel3D.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.u, 
                                _p.x, _p.y, _p.z, _p.s, _p.u, _p.neighbor);
        dudx = lsmpsVel3D.get_ddx();
        dudy = lsmpsVel3D.get_ddy();
        dudz = lsmpsVel3D.get_ddz();

        // Calculate the spatial differential for y-velocity 
        lsmpsVel3D.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.v, 
                                _p.x, _p.y, _p.z, _p.s, _p.v, _p.neighbor);
        dvdx = lsmpsVel3D.get_ddx();
        dvdy = lsmpsVel3D.get_ddy();
        dvdz = lsmpsVel3D.get_ddz();

        // Calculate the spatial differential for z-velocity 
        lsmpsVel3D.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.w, 
                                _p.x, _p.y, _p.z, _p.s, _p.w, _p.neighbor);
        dwdx = lsmpsVel3D.get_ddx();
        dwdy = lsmpsVel3D.get_ddy();
        dwdz = lsmpsVel3D.get_ddz();

		// Iteration for the velocity value
        for (int _i = 0; _i < _p.num; _i++){
			std::vector<double> velGrad(9);
            velGrad[0] = std::abs(dudx[_i]);
            velGrad[1] = std::abs(dudy[_i]);
            velGrad[2] = std::abs(dudz[_i]);

            velGrad[3] = std::abs(dvdx[_i]);
            velGrad[4] = std::abs(dvdy[_i]);
            velGrad[5] = std::abs(dvdz[_i]);
            
            velGrad[6] = std::abs(dwdx[_i]);
            velGrad[7] = std::abs(dwdy[_i]);
            velGrad[8] = std::abs(dwdz[_i]);

            // Find the maximum gradient value
            double maxGrad = *std::max_element(velGrad.begin(), velGrad.end());
			
			_lagStrNum[0] = maxGrad * Pars::dt;
			_lagStrNum[1] += _lagStrNum[0];
			_lagStrNum[2] = _lagStrNum[2] > _lagStrNum[0] ? _lagStrNum[2] : _lagStrNum[0];
		}
	#endif

    // Calculate the average value (put into the current evaluation container)
	_lagStrNum[0] = _lagStrNum[1] / _p.num;
    return;
}