/* ---------------- PROGRAM DESCRIPTION ----------------
	> Consisted of utilities function for main program
*/

#include "Utils.hpp"
#include "src/LSMPS/LSMPSa.hpp"

/**
 *  @brief  Start the simulation counter number. Set up the number of current 
 *  iteration digit.
 *         
 *  @param  step The iteration step.
*/
void simUtil::startCounter(int step){
	// Calculate the digit length of iteration step
	this->iterDigitLen = 1;
	if (step != 0){
		iterDigitLen = 1 + std::floor(std::log10(step));
	}
	return;
}

/**
 *  @brief  Printing the iteration number header to the console.
 *         
 *  @param  step The iteration step to be printed.
*/
void simUtil::printHeader(int step){
	// Calculate the digit length of iteration step
	int dig_len = 1;
	if (step != 0) dig_len = 1 + std::floor(std::log10(step));

	// Printing the iteration HEADER
	printf("\n\n+");
	for(int _i = 0; _i < 16 - dig_len/2; _i ++){
		printf("-");
	}
	printf(" Iteration Step %d ", step);
	for(int _i = 0; _i < 16 - dig_len/2 - dig_len%2; _i ++){
		printf("-");
	}
	printf("+\n");
	return;
}

/**
 *  @brief  Set the iteration label for file name.
 *         
 *  @param  step The iteration step.
 * 
 *  @return The iteration label name.
*/
std::string simUtil::saveName(const int step){
	std::string DataName;
	// Calculate the digit length of iteration step
	int dig_len = 1;
	if (step != 0) dig_len = 1 + std::floor(std::log10(step));
	
	// Add the leading zero
	int addDigit = Pars::max_dig_len - dig_len;
	for (int _spc = 0; _spc < addDigit; _spc++) DataName.append("0");

	// Combine the leading zero to step number
	DataName.append(std::to_string(step));
	return DataName;
}

/**
 *  @brief  Display the prediction time to finish the simulation.
 *         
 *  @param  step The current iteration step.
 *  @param  _currTime The current iteration computational time.
*/
void simUtil::predictCompTime(int step, double curr_comp_time){
	// The value of "curr_comp_time" is in second (s)
	// Internal variable
	int est_time_d, est_time_h, est_time_m; double est_time_s;
	est_time_s = curr_comp_time * double(Pars::max_iter - step - 1);
	
	// Calculate Day
	est_time_d = int(est_time_s / (24 * 60 * 60));
	est_time_s -= est_time_d * (24 * 60 * 60);
	// Calculate Hour
	est_time_h = int(est_time_s / (60 * 60));
	est_time_s -= est_time_h * (60 * 60);
	// Calculate Minute
	est_time_m = int(est_time_s / (60));
	est_time_s -= est_time_m * (60);
	
	// The simulation estimation is limited to 999 days of simulation
	printf("\n<!> Estimation time to finish run : %12.2f s", curr_comp_time*double(Pars::max_iter - step));
	if (est_time_d == 0){
		printf("\n<!> Estimation time to finish run :    %2dh %2dm %2ds", est_time_h, est_time_m, (int)est_time_s);
	}else{
		printf("\n<!> Estimation time to finish run :  %3ddays %2dhrs", est_time_d, est_time_h);
	}

	printf("\n");
	return;
}

/**
 *  @brief  Calculating the simulation residual and write into file.
 *         
 *  @param  _particle The particle data for residual calculation.
 *  @param  step The current iteration step.
*/
void simUtil::saveResidual(Particle &par, int step){
    // A prompt in stability evaluation
	printf("\nCalculating residuals ...\n");
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif
	
	// Internal variable
	std::vector<double> mass_cnsv(par.num,0.0);
	std::vector<double> vor(par.num,0.0);
    std::vector<double> vorx(par.num,0.0);
	std::vector<double> vory(par.num,0.0);
	std::vector<double> vorz(par.num,0.0);

	// Calculating the residual
	if (DIM == 2){
		// Internal tools
		// LSMPSa lsmpsa_du, lsmpsa_dv;
		LSMPSa lsmpsa_tool;

		// Calculate the differential du/dx
		lsmpsa_tool.set_LSMPS(par.x, par.y, par.s, par.u, par.neighbor);
		auto _dudx = lsmpsa_tool.get_ddx();
		auto _dudy = lsmpsa_tool.get_ddy();

		// Calculate the differential dv/dy
		lsmpsa_tool.set_LSMPS(par.x, par.y, par.s, par.v, par.neighbor);
		auto _dvdx = lsmpsa_tool.get_ddx();
		auto _dvdy = lsmpsa_tool.get_ddy();

		// Calculating vorticity and mass conservation
		for (int i = 0; i < par.num; i++){
			mass_cnsv[i] = _dudx[i] + _dvdy[i];
			vor[i] = par.vorticity[i] - (_dvdx[i] - _dudy[i]);
		}
	}
	else if (DIM == 3){
		// Internal tools
		// LSMPSa lsmpsa_du, lsmpsa_dv, lsmpsa_dw;
		LSMPSa lsmpsa_tool;

		// Calculate the differential du/dx
		lsmpsa_tool.set_LSMPS_3D(par.x, par.y, par.z, par.s, par.u, 
								 par.x, par.y, par.z, par.s, par.u, par.neighbor);
		auto _dudx = lsmpsa_tool.get_ddx();
		auto _dudy = lsmpsa_tool.get_ddy();
		auto _dudz = lsmpsa_tool.get_ddz();

		// Calculate the differential dv/dy
		lsmpsa_tool.set_LSMPS_3D(par.x, par.y, par.z, par.s, par.v, 
								 par.x, par.y, par.z, par.s, par.v, par.neighbor);
		auto _dvdx = lsmpsa_tool.get_ddx();
		auto _dvdy = lsmpsa_tool.get_ddy();
		auto _dvdz = lsmpsa_tool.get_ddz();

		// Calculate the differential dv/dy
		lsmpsa_tool.set_LSMPS_3D(par.x, par.y, par.z, par.s, par.w, 
								 par.x, par.y, par.z, par.s, par.w, par.neighbor);
		auto _dwdx = lsmpsa_tool.get_ddx();
		auto _dwdy = lsmpsa_tool.get_ddy();
		auto _dwdz = lsmpsa_tool.get_ddz();

		// Calculating vorticity

		for (int i = 0; i < par.num; i++){
			mass_cnsv[i] = _dudx[i] + _dvdy[i] + _dwdz[i];
			vorx[i] = par.vortx[i] - (_dwdy[i] - _dvdz[i]);
			vory[i] = par.vorty[i] - (_dudz[i] - _dwdx[i]);
			vorz[i] = par.vortz[i] - (_dvdx[i] - _dudy[i]);
		}
	}

	// Saving the residual field data
	if ((step % Pars::save_inv) == 0){
		// Initialize the writter
		std::ofstream data;
		std::string name;
		name  = "output/residual_";
		simUtil util_step;
		std::string number = util_step.saveName(step);
		name += number + ".csv";
		
		// Write data
		data.open(name);
		
		if (DIM == 2){
			// Write header
			data << "xp,yp,sp,vor,mass_cnsv\n";
			// Write data content
			for (int i = 0; i < par.num; i++){
			data << ""  << par.x[i]
				<< "," << par.y[i]
				<< "," << par.s[i]
				<< "," << vor[i]
				<< "," << mass_cnsv[i]
				<< "\n";
			}
		}
		else if (DIM == 3){
			// Write header
			data << "xp,yp,zp,sp,vorx,vory,vorz,mass_cnsv\n";
			// Write data content
			for (int i = 0; i < par.num; i++){
			data << ""  << par.x[i]
				<< "," << par.y[i]
				<< "," << par.z[i]
				<< "," << par.s[i]
				<< "," << vorx[i]
				<< "," << vory[i]
				<< "," << vorz[i]
				<< "," << mass_cnsv[i]
				<< "\n";
			}
		}
		data.close();
	}


	// Calculate the global residual
	double resMC = 0.0;		// Residual of mass conservation
	double resVor = 0.0;	// Residual of vorticity
	
	// Residual Calculation
	#if (DIM == 2)
		for (int i = 0; i < par.num; i++){
			// Residual of mass conservation E_mc = sum(abs(div(U))*Area) / (U_inf * D)
			resMC += std::abs(mass_cnsv[i]) * par.s[i]*par.s[i];

			// Residual of vorticity E_vor = sum(|curl(U) - omega|^2 * Area)  / (U_inf^2)
			resVor += vor[i]*vor[i] * par.s[i]*par.s[i];
		}
		
		// Normalize the value
		resMC /= (Pars::U_inf * Pars::Df);
		resVor /= (Pars::U_inf * Pars::U_inf);
	
	#elif (DIM == 3)
		for (int i = 0; i < par.num; i++){
			// Residual of mass conservation E_mc = sum(abs(div(U))*Vol) / (U_inf * D^2)
			resMC += std::abs(mass_cnsv[i]) * par.s[i]*par.s[i]*par.s[i];

			// Residual of vorticity E_vor = sum(|curl(U) - omega|^2 * Vol)  / (U_inf^2*D)
			resVor += (vorx[i]*vorx[i] + vory[i]*vory[i] + vorz[i]*vorz[i]) * par.s[i]*par.s[i]*par.s[i];
		}
		
		// Normalize the value
		resMC /= (Pars::U_inf * Pars::Df * Pars::Df);
		resVor /= (Pars::U_inf * Pars::U_inf * Pars::Df);
	#endif
	
	// Saving the residual history data
	std::ofstream _write;
	std::string fileName = "output/Residual_History";
	if      (Pars::opt_start_state == 0) fileName += ".csv";
	else if (Pars::opt_start_state == 1) fileName += "_Resume.csv";
	
	// Determine the current step is initialization or not
	bool init = false;
	if      (Pars::opt_start_state == 0) init = (step == 0);
	else if (Pars::opt_start_state == 1) init = (step == (Pars::resume_step + 1));

	// Save the header
	if (init == true){
		_write.open(fileName);
		_write << "time,massRes,vorticityRes\n";
		_write.close();
	}
	
	// Save the data
	_write.open(fileName, std::ofstream::out | std::ofstream::app);
	_write <<  "" << step*Pars::dt 
			<< "," << resMC
			<< "," << resVor
			<< "\n";
	_write.close();


	// Particle generation initialization summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Stability evaluation comp. time.   [%f s]\n", _time);
    
    return;
}

void simUtil::addVorPertubation(Particle &p){
	// Add the peturbation right in front of the object
	// with value of vorticity = |omega|_max * 10 ^-3

	// The peturbation function is a section of cosine function
	//    omega_ptb = (K/2)*(1 + cos(pi*(r/R))) ; (r < R)
	//  Where : 
	//    > K : vorticity peak multiplication
	//    > Sc: vorticity scaling
	//    > R : peturbation radius -> (sigma * Rs)
	//    > r : radius from peturbation location
	//    > sigma : local max particle size
	//    > Rs : peturbation radius scale

	// Set the peturbation constant
	double K, Sc, R, Rs;
	Rs = Pars::r_sup;
	Sc = 0.05;

	// Set the peturbation location
	double ptbCoord[DIM];
	ptbCoord[0] = -0.7;    // x location
	ptbCoord[1] = -0.6;    // y location
	#if DIM == 3
	ptbCoord[2] = 0.0;     // z location
	#endif

	// Evaluate the box size
	double evalBoxSize = Pars::sigma * 16;
	double locMaxParSize = Pars::sigma;		// The local maximum particle size
	double vorMaxVal = 0.0;		// The local maximum particle size

	// Find the local maximum particle size and maximum vorticity
	for (int i = 0; i < p.num; i++){
		// Update the maximum value of vorticity
		if (std::abs(p.vorticity[i]) > vorMaxVal) vorMaxVal = std::abs(p.vorticity[i]);

		// Check whether the particle is outside evaluation box
		if (p.x[i] > ptbCoord[0] + evalBoxSize ||
			p.x[i] < ptbCoord[0] - evalBoxSize)
		continue;

		if (p.y[i] > ptbCoord[1] + evalBoxSize ||
			p.y[i] < ptbCoord[1] - evalBoxSize)
		continue;

		#if DIM == 3
		if (p.z[i] > ptbCoord[2] + evalBoxSize ||
			p.z[i] < ptbCoord[2] - evalBoxSize)
		continue;
		#endif

		// Evaluate the particle size
		if (p.s[i] > locMaxParSize) locMaxParSize = p.s[i];
	}

	// Update the evaluation box size and peturbation parameter
	evalBoxSize = locMaxParSize * Rs * 1.2;
	R = locMaxParSize * Rs;
	K = vorMaxVal * Sc;

	// Calculate Set the peturbation value
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		// [CHECK 1] Check whether the particle is outside evaluation box
		if (p.x[i] > ptbCoord[0] + evalBoxSize ||
			p.x[i] < ptbCoord[0] - evalBoxSize)
		continue;

		if (p.y[i] > ptbCoord[1] + evalBoxSize ||
			p.y[i] < ptbCoord[1] - evalBoxSize)
		continue;

		#if DIM == 3
		if (p.z[i] > ptbCoord[2] + evalBoxSize ||
			p.z[i] < ptbCoord[2] - evalBoxSize)
		continue;
		#endif

		// Update the vorticity
		double r2 = 0.0;
		double dx, dy, dz;
		// Add the x component
			dx = p.x[i] - ptbCoord[0];
			r2 += dx*dx;
		// Add the y component
			dy = p.y[i] - ptbCoord[1];
			r2 += dy*dy;
		// Add the z component
		#if DIM == 3
			dz = p.z[i] - ptbCoord[2];
			r2 += dz*dz;
		#endif
		double r = std::sqrt(r2);

		// [CHECK 2] Check the location of particle related to the peturbation location
		if (r > R) continue;

		// Calculate the local vorticity
		#if (DIM == 2)
			// Peturbation for 2D space (update on vort z)
			double ptbVor = K/2.0 * (1 + std::cos(M_PI * r / R));
			p.vorticity[i] += ptbVor;
			p.gz[i] = p.vorticity[i] * p.s[i]*p.s[i];
		#elif (DIM == 3)
			// Peturbation for 3D space (update on vort z)
			double rTh = std::sqrt(dx*dx + dy*dy);
			double ptbVorTh = 0.5 * (1 + std::cos(M_PI * rTh / std::sqrt(R*R - dz*dz)));
			double ptbVorZ = 0.5 * (1 + std::cos(M_PI * dz / R));
			p.vortz[i] += K*ptbVorZ*ptbVorTh;
		#endif
	}

	return;
}
void simUtil::calculate_save_enstrophy(Particle &par, int it){
	// Initialized internal variable
	double E = 0.0;
	double G = 0.0;
	std::ofstream writer;
	std::string fileName = "output/enstrophy_history.csv";

	// Calculate the enstrophy integration
	// #pragma omp parallel for
	for (int i = 0; i < par.num; i++){
		double vol = std::pow(par.s[i], DIM);
		#if DIM == 2
			double _E = (par.vorticity[i]*par.vorticity[i]*vol);
			// double _G = (par.vorticity[i]*vol);
			// G += _G;
		#elif DIM == 3
			double _E = (par.vortx[i]*par.vortx[i] + par.vorty[i]*par.vorty[i] + par.vortz[i]*par.vortz[i]) * vol;
		#endif
		E += _E;
	}

	#if DIM == 2
		// Calculate the error of Enstrophy
		double T = it*Pars::dt;
		double L = Pars::L_LO;
		double Re = Pars::RE;
		double w0 = Re*Pars::NU / (M_PI*L*L);
		double f_t = Re / (4*M_PI*T*w0 + Re);
		double E_a = (M_PI/2.0)*(w0*w0)*(L*L)*f_t;
	#endif

	// Additional of global maximum vorticity and minimum vorticity
	double vorMax = *(std::max_element(par.vorticity.begin(),par.vorticity.end()));
	double vorMin = *(std::min_element(par.vorticity.begin(),par.vorticity.end()));

	// // Write the data header
	// if (it == 0){
	// 	writer.open(fileName);
	// 	writer << "time,vorMin,vorMax,EnstA,Enst,errE,Circ,errC\n";
	// 	writer.close();
		
	// 	// Update the initial enstrophy
	// 	this->E_0 = E;
	// 	// Update the initial circulation
	// 	this->G_0 = G;
	// }

	// // Write the element data
	// writer.open(fileName, std::ofstream::out | std::ofstream::app);
	// writer << it * Pars::dt
	// 	   << "," << vorMin
	// 	   << "," << vorMax
	// 	   << "," << E_a
	// 	   << "," << E
	// 	   << "," << (E - E_a)/E_a
	// 	   << "," << G
	// 	   << "," << (G - this->G_0)/this->G_0
	// 	   << "\n";
	// writer.close();
	
	// Write the data header
	if (it == 0){
		writer.open(fileName);
		writer << "time,vorMin,vorMax,Enst,errE\n";
		writer.close();
		
		// Update the initial enstrophy
		this->E_0 = E;
		// // Update the initial circulation
		// this->G_0 = G;
	}

	// Write the element data
	writer.open(fileName, std::ofstream::out | std::ofstream::app);
	writer << it * Pars::dt
		   << "," << vorMin
		   << "," << vorMax
		   << "," << E
		   << "," << (E - this->E_0)/this->E_0
		//    << "," << G
		//    << "," << (G - this->G_0)/this->G_0
		   << "\n";
	writer.close();
	return;
}

void simUtil::calculate_save_dissipation(Particle &par, int it){
	// Initialized internal variable
	double J = 0.0;
	double G = 0.0;
	std::ofstream writer;
	std::string fileName = "output/disipation_history.csv";

	// Here we will calculate the total circulation and also the total dissipation

	// Calculate the total circulation and angular momentum
	// #pragma omp parallel for
	for (int i = 0; i < par.num; i++){
		double vol = std::pow(par.s[i], DIM);
		double _G = (par.vorticity[i]*vol);
		double _J = (par.vorticity[i]*(par.x[i]*par.x[i] + par.y[i]*par.y[i])*vol);
		J += _J;
		G += _G;
	}
	
	// Write the data header
	if (it == 0){
		writer.open(fileName);
		writer << "time,gTot,dG,jTot,dJ,ReV\n";
		writer.close();

		//  Note :
		// gTot = G
		// dG = G - G_0
		// aTot = J
		// dJ = J - J_0
		// ReV = G^2 * t/ (dJ)
		
		// Update the initial enstrophy
		this->J_0 = J;
		this->G_0 = G;
	}

	// Write the element data
	writer.open(fileName, std::ofstream::out | std::ofstream::app);
	writer << it * Pars::dt
		   << "," << G
		   << "," << (G - this->G_0)
		   << "," << J
		   << "," << (J - this->J_0)
		   << "," << 4*G*G*it*Pars::dt / (J - this->J_0)
		   << "\n";
	writer.close();

	return;
}

void simUtil::calculate_save_kinetic(Particle &par, int it){
	// Initialized internal variable
	double E = 0.0;
	std::ofstream writer;
	std::string fileName = "output/kinetic_history.csv";

	// Calculate the enstrophy integration
	// #pragma omp parallel for
	for (int i = 0; i < par.num; i++){
		double U2 = par.u[i]*par.u[i] + par.v[i]*par.v[i];
		#if DIM == 3
			U2 += par.w[i]*par.w[i];
		#endif
		double vol = std::pow(par.s[i], DIM);
		double _E = 0.5 * U2 * vol;
		E += _E;
	}
	
	// Write the data header
	if (it == 0){
		writer.open(fileName);
		writer << "time,kinetic\n";
		writer.close();
	}

	// Write the element data
	writer.open(fileName, std::ofstream::out | std::ofstream::app);
	writer << it * Pars::dt
		   << "," << E
		   << "\n";
	writer.close();

	return;
}

// Note! Still limited for 2D simulation calculation
// The calculation of L-Norm Error is given below

/**
 * Supposed that a function f have 
 * numerical solution of fn and analytical solution of fa
 *  > The L2 Norm Error (Global representation)
 *                /  sum_(i)^(N) {(fn_i - fa_i)^2}  \
 *      L2 = sqrt | ------------------------------- |
 *                \      sum_(i)^(N) {fa_i^2}       /
 * 
 *  > The L_inf Norm Error (Maximum Value Representation)
 *                   /  max {fn_i - fa_i}  \
 *      L_inf =      | ------------------- |
 *                   \     max {fa_i}      /
 *      
*/
double simUtil::calculate_L2_Norm(const Particle &par, const std::vector<double> &f_n, const std::vector<double> &f_a){
	// Initialize the norm
	double SSD = 0.0e0;		// The sum of square different between numerical and analytical
	double SAS = 0.0e0;		// The sum of analytical value square

	// The domain confining
	const double cut_off_dist = Pars::sigma * 5;
	const double lowerX = -Pars::xdom + cut_off_dist;
	const double upperX =  Pars::lxdom - Pars::xdom - cut_off_dist;
	const double lowerY = -0.5*Pars::lydom + cut_off_dist;
	const double upperY =  0.5*Pars::lydom - cut_off_dist;

	// Calculate the norm error
	for (int i = 0; i < par.num; i++){
		// Only evaluate particle inside the bound
		if (par.x[i] > lowerX && par.x[i] < upperX && 
			par.y[i] > lowerY && par.y[i] < upperY )
		{
			double diff = f_n[i] - f_a[i];
			
			// Calculate L2 Norm evaluation
			SSD += diff*diff;
			SAS += f_a[i]*f_a[i];
		}
	}
	
	// Calculate the finals
	double L2N = std::sqrt(SSD/SAS);

	return L2N;
}

void simUtil::calculate_L_Norm(std::vector<double>& Norm, const Particle &par, const std::vector<double> &f_n, const std::vector<double> &f_a){
	// Initialize the norm
	double SSD = 0.0e0;		// The sum of square different between numerical and analytical
	double SAS = 0.0e0;		// The sum of analytical value square

	// Here also calculates the maximum error
	double SSD_max = 0.0e0;		// The max value of square different between numerical and analytical
	double SAS_max = 0.0e0;		// The max value of analytical value square

	// The domain confining
	const double cut_off_dist = Pars::sigma * 5;
	const double lowerX = -Pars::xdom + cut_off_dist;
	const double upperX =  Pars::lxdom - Pars::xdom - cut_off_dist;
	const double lowerY = -0.5*Pars::lydom + cut_off_dist;
	const double upperY =  0.5*Pars::lydom - cut_off_dist;

	// Calculate the norm error
	for (int i = 0; i < par.num; i++){
		// Only evaluate particle inside the bound
		if (par.x[i] > lowerX && par.x[i] < upperX && 
			par.y[i] > lowerY && par.y[i] < upperY )
		{
			double diff = f_n[i] - f_a[i];
			
			// Calculate L2 Norm evaluation
			SSD += diff*diff;
			SAS += f_a[i]*f_a[i];

			// Calculate L_inf Norm evaluation
			SSD_max = std::max(SSD_max, std::abs(diff));
			SAS_max = std::max(SAS_max, std::abs(f_a[i]));
		}
	}
	
	// Calculate the finals
	double L2N = std::sqrt(SSD/SAS);
	double LinfN = SSD_max/SAS_max;

	// Return the value by reference
	Norm.clear();
	Norm.push_back(L2N);
	Norm.push_back(LinfN);

	return;
}

#include "src/initialization/initialization.hpp"
/**
 *  @brief  Calculate the L2 Norm error of velocity. A global error calculation evaluation
 *         
 *  @param  par  The particle data container to be evaluated.
 *  @param  it   The iteration step on the current evaluation.
 *  @param  _func Type of function to be compared.
 *                1:= Perlman vorticity,
 *                2:= Lamb Oseen vorticity
 *  @param  _type Type of evaluation and data saving
 *                0:= Evaluate L2 Norm using all component
 *                1:= Evaluate L2 Norm using total velocity 
 *                2:= Evaluate each velocity component
 *  @return No return.
*/
void simUtil::velocity_L2Norm(Particle &par, int it, int _func, int _type){
	// Initialized internal variable
	std::ofstream writer;
	std::string fileName = "output/velocity_L2_history.csv";

	// Calculate the analytical velocity data
	initialization init_tool;
	if (_func == 1) init_tool.perlman_velocity_solution(par, 1);
	else if (_func == 2) init_tool.lamb_oseen_solution(par, it*Pars::dt, 0);

	// The value of L2 norm error of each value
	double L2N_Vel;
	double L2N_u;
	double L2N_v;
	double L2N_w;

	// The value of Linf norm error of each value
	double Linf_N_Vel;
	double Linf_N_u;
	double Linf_N_v;
	double Linf_N_w;

	// Evaluation for total vorticity
	// ******************************
	if (_type == 0 || _type == 1){
		// Calculate the L2 norm
		double SSD = 0.0e0;		// The sum of square different between numerical and analytical
		double SAS = 0.0e0;		// The sum of analytical value square

		// Calculate the Linf norm
		double SSD_max = 0.0e0;		// The max of square different between numerical and analytical
		double SAS_max = 0.0e0;		// The max of analytical value square

		// The domain confining
		const double cut_off_dist = Pars::sigma * 5;
		const double lowerX = -Pars::xdom + cut_off_dist;
		const double upperX =  Pars::lxdom - Pars::xdom - cut_off_dist;
		const double lowerY = -0.5*Pars::lydom + cut_off_dist;
		const double upperY =  0.5*Pars::lydom - cut_off_dist;

		// Calculate the norm error
		for (int i = 0; i < par.num; i++){
			// Only evaluate particle inside the bound
			if (par.x[i] > lowerX && par.x[i] < upperX && 
				par.y[i] > lowerY && par.y[i] < upperY )
			{
				// Calculate the velocity vector difference
				double du = par.u_a[i] - par.u[i];
				double dv = par.v_a[i] - par.v[i];

				// Calculate the L2 Norm component
				double diff = du*du + dv*dv;
				double analytic = par.u_a[i]*par.u_a[i] + par.v_a[i]*par.v_a[i];

				// Additional calculation for 3D (Not in calculation)
				#if DIM == 3
					double dw = par.w_a[i] - par.w[i];
					diff += dw*dw;
					analytic += par.w_a[i]*par.w_a[i];
				#endif

				// Calculate L2 Norm evaluation
				SSD += diff;
				SAS += analytic;

				// Calculate L_inf Norm evaluation
				SSD_max = std::max(SSD_max, sqrt(diff));
				SAS_max = std::max(SAS_max, sqrt(analytic));
			}
		}
		
		// Calculate the finals
		L2N_Vel = std::sqrt(SSD/SAS);
		Linf_N_Vel = SSD_max/SAS_max;
	}
	if (_type == 0 || _type == 2){
		// Compact value of L norm error of each value
		std::vector<double> L_N_Vel;
		std::vector<double> L_N_u;
		std::vector<double> L_N_v;
		std::vector<double> L_N_w;

		// Calculate the L2 norm
		// L2N_u = this->calculate_L2_Norm(par, par.u, par.u_a);
		// L2N_v = this->calculate_L2_Norm(par, par.v, par.v_a);
		this->calculate_L_Norm(L_N_u, par, par.u, par.u_a);
		this->calculate_L_Norm(L_N_v, par, par.v, par.v_a);
		
		// Assign the value
		L2N_u = L_N_u[0];
		L2N_v = L_N_v[0];
		Linf_N_u = L_N_u[1];
		Linf_N_v = L_N_v[1];

		// Additional calculation for 3D (Not in calculation)
		#if DIM == 3
			this->calculate_L_Norm(L_N_w, par, par.w, par.w_a);
			L2N_w = L_N_w[0];
			Linf_N_w = L_N_w[1];
		#endif
	}

	
	// Write the data
	// **************
	if (_type == 0){
		// Write the data header
		if (it == 0){
			writer.open(fileName);
			writer << "time,L2N_u,L2N_v,L2N_Vel,LIN_u,LIN_v,LIN_Vel\n";
			writer.close();
		}

		// Write the element data
		writer.open(fileName, std::ofstream::out | std::ofstream::app);
		writer << it * Pars::dt
			<< "," << L2N_u
			<< "," << L2N_v
			<< "," << L2N_Vel
			<< "," << Linf_N_u
			<< "," << Linf_N_v
			<< "," << Linf_N_Vel
			<< "\n";
		writer.close();
	}

	else if (_type == 1){
		// Write the data header
		if (it == 0){
			writer.open(fileName);
			writer << "time,L2N_Vel,LIN_Vel\n";
			writer.close();
		}

		// Write the element data
		writer.open(fileName, std::ofstream::out | std::ofstream::app);
		writer << it * Pars::dt
			<< "," << L2N_Vel
			<< "," << Linf_N_Vel
			<< "\n";
		writer.close();
	}

	else if(_type == 2){
		// Write the data header
		if (it == 0){
			writer.open(fileName);
			writer << "time,L2N_u,L2N_v,LIN_u,LIN_v\n";
			writer.close();
		}

		// Write the element data
		writer.open(fileName, std::ofstream::out | std::ofstream::app);
		writer << it * Pars::dt
			<< "," << L2N_u
			<< "," << L2N_v
			<< "," << Linf_N_u
			<< "," << Linf_N_v
			<< "\n";
		writer.close();
	}

	return;
}


/**
 *  @brief  Calculate the L2 Norm error of velocity. A global error calculation evaluation
 *         
 *  @param  par  The particle data container to be evaluated.
 *  @param  it   The iteration step on the current evaluation.
 *  @param  _func Type of function to be compared.
 *                1:= Perlman vorticity,
 *                2:= Lamb Oseen vorticity
 *  @return No return.
*/
void simUtil::vorticity_L2Norm(Particle &par, int it, int _func){
	// The function is already calculated !!!
	// Initialized internal variable
	std::ofstream writer;
	std::string fileName = "output/vorticity_L2_history.csv";
	
	// Calculate the finals	
	std::vector<double> LNorm;
	// double L2N = this->calculate_L2_Norm(par, par.vorticity, par.vort_a);
	this->calculate_L_Norm(LNorm, par, par.vorticity, par.vort_a);
	double L2N = LNorm[0];
	double Linf_N = LNorm[1];
	
	// Write the data header
	if (it == 0){
		writer.open(fileName);
		writer << "time,L2N_Vort,LinfN_Vort\n";
		writer.close();
	}

	// Write the element data
	writer.open(fileName, std::ofstream::out | std::ofstream::app);
	writer << it * Pars::dt
		<< "," << L2N
		<< "," << Linf_N
		<< "\n";
	writer.close();

	return;
}

/**
 *  @brief A utility function to take only the particle with nonzero vorticity
 *  and the particle close to it.
 *  PS: Currently only collect data
 *      1. Computational Parameter (size, neighbor, active flag)
 *      2. Coordinate
 *      3. Physical Properties (velocity, vorticity)
 *      4. Grid Node parameter [v]
 *      5. Body parameter [x] -> Still developed
 *  
 *  @param _particle The particle data for trim domain evaluation.
 */
void simUtil::trimParticleDomain(Particle &particle){
    // [!] Note that the nonzero vorticity particle must already have the active sign !!!
	
	// Diffusion calculation computational time manager
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // Internal variables
	Particle& p = particle;      			// The data of all particle in the domain
	std::vector<int> _index;                // The index list of particle inside 'active particle + buffer zone' (active -> original)
	std::unordered_map<int,int> _inv_index; // The index container to take the active particle index (original -> active)
	std::vector<bool> _flag(particle.num, false);      // The flag list of evaluated particle, At initial still no evaluated particle list

	// PROCEDURE 1:
	// ************
	// Collecting the only active particle data
	// Assign the index of active particles
	for (int ID = 0; ID < p.num; ID++)
	{
		// Check if the particle ID is active
		if (p.isActive[ID] == true)
		{
			// [1] The current ID particle
			if(!Pars::flag_ngh_include_self && _flag[ID] == false)
			{
				// Assign ID into index list if not in the list (if flag == false)
				_inv_index[ID] = _index.size();
				_index.push_back(ID);    // Assign ID into the index list
				_flag[ID] = true;        // Update the evaluation flag of particle ID
			}

			// [2] The neighbor particles of current ID particle
			for (auto _ngh_ID : p.neighbor[ID])
			{
				if(_flag[_ngh_ID] == false)
				{
					// Assign into index list if not in the list (if flag == false)
					_inv_index[_ngh_ID] = _index.size();
					_index.push_back(_ngh_ID);  // Assign _ngh_ID particle into the index list
					_flag[_ngh_ID] = true;      // Update the evaluation flag of particle _ngh_ID
				}
			}
		}
	}

	// Generate the new particle
	Particle _particle;

	// Particle resize
	// ***************
	_particle.num = _index.size();
	_particle.s.clear(); _particle.s.resize(_particle.num);
	_particle.isActive.clear(); _particle.isActive.resize(_particle.num);
	_particle.neighbor.clear(); _particle.neighbor.resize(_particle.num);
	if (Pars::opt_init_particle == 5){
		_particle.level.clear(); _particle.level.resize(_particle.num);
		_particle.nodeID.clear(); _particle.nodeID.resize(_particle.num);
	}
	
	// Coordinate
	// **********
	_particle.x.clear(); _particle.x.resize(_particle.num);
	_particle.y.clear(); _particle.y.resize(_particle.num);
	#if DIM == 3
		_particle.z.clear(); _particle.z.resize(_particle.num);
	#endif

	// Properties
	// **********
	// Vorticities
	#if DIM == 3
		_particle.vortx.clear(); _particle.vortx.resize(_particle.num);
		_particle.vorty.clear(); _particle.vorty.resize(_particle.num);
		_particle.vortz.clear(); _particle.vortz.resize(_particle.num);
	#endif
	_particle.vorticity.clear(); _particle.vorticity.resize(_particle.num);

	// Velocities
	_particle.u.clear(); _particle.u.resize(_particle.num);
	_particle.v.clear(); _particle.v.resize(_particle.num);
	#if DIM == 3
		_particle.w.clear(); _particle.w.resize(_particle.num);
	#endif

	// Store the particle data of each particle inside _index list
	#pragma omp parallel for
	for (int i = 0; i < _particle.num; i++)
	{
		const int &_ID = _index[i];
		_particle.s[i] = (p.s[_ID]);
		_particle.isActive[i] = (p.isActive[_ID]);

		_particle.x[i] = (p.x[_ID]);
		_particle.y[i] = (p.y[_ID]);
		_particle.u[i] = (p.u[_ID]);
		_particle.v[i] = (p.v[_ID]);
		_particle.vorticity[i] = (p.vorticity[_ID]);

		#if DIM == 3
		_particle.z[i] = (p.z[_ID]);
		_particle.w[i] = (p.w[_ID]);
		_particle.vortx[i] = (p.vortx[_ID]);
		_particle.vorty[i] = (p.vorty[_ID]);
		_particle.vortz[i] = (p.vortz[_ID]);
		#endif
		
		for (int oriNghID : p.neighbor[_ID]){
			// Check whether the particle at neighbor ID of the original data is existed in trimmed particle data
			if (_inv_index.count(oriNghID) != 0){
				_particle.neighbor[i].push_back(_inv_index.at(oriNghID));
			}
		}

	}
	
	// Update the node ID and level catagory only for the node grid based particle initialization
	if (Pars::opt_init_particle == 5){
		#pragma omp parallel for
		for (int i = 0; i < _particle.num; i++)
		{
			const int &_ID = _index[i];
			_particle.level[i] = (p.level[_ID]);
			_particle.nodeID[i] = (p.nodeID[_ID]);
		}
	}

	// Update the particle
	particle = _particle;

    // *Note (3D): At this point dont need to calculate the vorticity absolute value,
    //         because it will be calculated in the next section (particle redistribution)

    // Display local computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Collecting active particle:        [%f s]\n", _time);

    return;
}

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// The function list for LSMPS CALCULATION
void Util_LSMPS_gaussian_1(Particle& par){
	// Prompt
	std::cout << "Start evaluating the gaussian!\n";

	// Parameter
	double c = 10.0e0;		// The function scaling
	double a = 1.0e0;		// The width factor

	// Resize the data
	par.F_a.resize(par.num);
	par.Fx_a.resize(par.num);
	par.Fy_a.resize(par.num);
	par.Fx2_a.resize(par.num);
	par.Fy2_a.resize(par.num);

	// Calculate the function component
	#pragma omp parallel for
	for (int i = 0; i < par.num; i++){
		// A reference
		double& X = par.x[i];
		double& Y = par.y[i];

		// Calculation Component
		double F = exp(-a*(X*X + Y*Y));
		double F_d = -2*a*F;

		// Calculation Final Result
		par.F_a[i]   = c*( F );
		par.Fx_a[i]  = c*( X*F_d );
		par.Fy_a[i]  = c*( Y*F_d );
		par.Fx2_a[i] = c*( F_d * (1 - 2*a*X*X) );
		par.Fy2_a[i] = c*( F_d * (1 - 2*a*Y*Y) );
	}

	std::cout << "[>] Done current step!\n";

	return;
}

void Util_LSMPS_gaussian_2(Particle& par){
	// Prompt
	std::cout << "Start evaluating the gaussian!\n";

	// Parameter
	double c = 1.0e0;		// The function scaling
	double a = 5.0e0;		// The width factor
	double R = 2.0;			// The radius of calculation
	const double k = log(2.0) / 4.0;	// A factor

	// // Intermediate parameter
	// std::vector<double> F(par.num), G(par.num), H(par.num);
	// std::vector<double> F_d(par.num), G_d(par.num), H_d(par.num);

	// Resize the data
	par.F_a.resize(par.num);
	par.Fx_a.resize(par.num);
	par.Fy_a.resize(par.num);
	par.Fx2_a.resize(par.num);
	par.Fy2_a.resize(par.num);

	// Calculate the function component
	#pragma omp parallel for
	for (int i = 0; i < par.num; i++){
		// A reference
		double& X = par.x[i];
		double& Y = par.y[i];
		double sqrt_x2y2 = sqrt(X*X + Y*Y);

		// Calculation Component
		double F = exp(-a*(sqrt_x2y2 - R)*(sqrt_x2y2 - R));
		double F_d = -2*a*F*(1 - R/sqrt_x2y2);
		double G = exp(-k*(X*X + Y*Y));
		double G_d = -2*k*G;
		double H = a + k - a*R/(sqrt_x2y2);
		double H_d = a*R/(sqrt_x2y2*sqrt_x2y2*sqrt_x2y2);

		// Calculation Final Result
		par.F_a[i]   = c*( Y*F*G );
		par.Fx_a[i]  = c*( -2*X*Y*F*G*H );
		par.Fy_a[i]  = c*( F*G - 2*Y*Y*F*G*H );
		par.Fx2_a[i] = c*( -2*Y*F*G*H - 2*X*X*Y*(F_d*G*H + F*G_d*H + F*G*H_d) );
		par.Fy2_a[i] = c*( -6*Y*F*G*H - 2*Y*Y*Y*(F_d*G*H + F*G_d*H + F*G*H_d) );
	}

	std::cout << "[>] Done current step!\n";

	return;
}

void Util_LSMPS_gaussian_3(Particle& par){
	// Prompt
	std::cout << "Start evaluating the gaussian!\n";

	// Parameter
	double c = 5.0e0;		// The function scaling

	// Resize the data
	par.F_a.resize(par.num);
	par.Fx_a.resize(par.num);
	par.Fy_a.resize(par.num);
	par.Fx2_a.resize(par.num);
	par.Fy2_a.resize(par.num);

	// Calculate the function component
	#pragma omp parallel for
	for (int i = 0; i < par.num; i++){
		// A reference
		double& X = par.x[i];
		double& Y = par.y[i];

		// Calculation Component
		double F = exp(-(X*X + Y*Y));
		double F_d = -2*F;

		// Calculation Final Result
		par.F_a[i]   = c*( X*Y*F );
		par.Fx_a[i]  = c*( Y*F*(1-2*X*X) );
		par.Fy_a[i]  = c*( X*F*(1-2*Y*Y) );
		par.Fx2_a[i] = c*( (4*X*X - 6)*X*Y*F );
		par.Fy2_a[i] = c*( (4*Y*Y - 6)*X*Y*F );
	}

	std::cout << "[>] Done current step!\n";

	return;
}

void Util_save_LSMPS_data(Particle & par, std::string name){
	// Basic parameter
	std::ofstream write;
	std::string name1;
	std::cout << "Start saving the data LSMPS\n";

	// // Print header
	// name1 = "output/Analytical_" + name + ".csv";
	// write.open(name1.c_str());
	// write <<  "" << "x"
	//       << "," << "y"
	// 	  << "," << "s"
	// 	  << "," << "f"
	// 	  << "," << "fx"
	// 	  << "," << "fy"
	// 	  << "," << "fx2"
	// 	  << "," << "fy2"
	// 	  << "," << "flap"
	// 	  << "\n" ;

	// // Print data
	// for (int i = 0; i < par.num; i++){
	// 	write <<  "" << par.x[i]
	// 		  << "," << par.y[i]
	// 		  << "," << par.s[i]
	// 		  << "," << par.F_a[i]
	// 		  << "," << par.Fx_a[i]
	// 		  << "," << par.Fy_a[i]
	// 		  << "," << par.Fx2_a[i]
	// 		  << "," << par.Fy2_a[i]
	// 		  << "," << par.Fx2_a[i] + par.Fy2_a[i]
	// 		  << "\n" ;
	// }

	// write.close();

	name1 = "output/LSMPS_" + name + ".csv";
	write.open(name1.c_str());
	write <<  "" << "x"
	      << "," << "y"
	      << "," << "s"
		  << "," << "f"
		  << "," << "fx"
		  << "," << "fy"
		  << "," << "fx2"
		  << "," << "fy2"
		  << "," << "flap"
		  << "\n" ;

	// Print data
	for (int i = 0; i < par.num; i++){
		write <<  "" << par.x[i]
			  << "," << par.y[i]
			  << "," << par.s[i]
			  << "," << par.F_a[i]
			  << "," << par.Fx[i]
			  << "," << par.Fy[i]
			  << "," << par.Fx2[i]
			  << "," << par.Fy2[i]
			  << "," << par.Fx2[i] + par.Fy2[i]
			  << "\n" ;
	}

	// write.close();
	
	std::cout << "[>] Done saving step\n";
	return;
}

void Util_save_LSMPS(Particle & par, std::string name){
	// Basic parameter
	std::ofstream write;
	std::string name1;

	// Calculation of Norm, Put into a dat document
	simUtil util;
	std::vector<double> tempNorm; // L2N, LInf
	// double L2N_f,   LInf_f;
	double L2N_fx,  LInf_fx;
	double L2N_fy,  LInf_fy;
	double L2N_fx2, LInf_fx2;
	double L2N_fy2, LInf_fy2;
	double L2N_lap, LInf_lap;

	// Calculate the norm error
	// util.calculate_L_Norm(tempNorm, par, par.F, par.F_a);
	util.calculate_L_Norm(tempNorm, par, par.Fx, par.Fx_a); 
	L2N_fx = tempNorm[0]; LInf_fx = tempNorm[1];
	
	util.calculate_L_Norm(tempNorm, par, par.Fy, par.Fy_a);
	L2N_fy = tempNorm[0]; LInf_fy = tempNorm[1];

	util.calculate_L_Norm(tempNorm, par, par.Fx2, par.Fx2_a);
	L2N_fx2 = tempNorm[0]; LInf_fx2 = tempNorm[1];

	util.calculate_L_Norm(tempNorm, par, par.Fy2, par.Fy2_a);
	L2N_fy2 = tempNorm[0]; LInf_fy2 = tempNorm[1];

	// Laplacian
	std::vector<double> lap(par.num), lap_a(par.num);
	for (int i = 0; i < par.num; i++){
		lap[i] = par.Fx2[i] + par.Fy2[i];
		lap_a[i] = par.Fx2_a[i] + par.Fy2_a[i];
	}
	util.calculate_L_Norm(tempNorm, par, lap, lap_a);
	L2N_lap = tempNorm[0]; LInf_lap = tempNorm[1];

	std::cout << "[>] Done calculate norm\n";

	// Print Norm
	name1 = "output/" + name + ".dat";
	write.open(name1.c_str());
	write << "The L norm summary\n";
	write << " >> DATA: "+ name + "\n";
	write << "Particle : " << par.num << "\n";

	write << "\nThe L2 Norm\n"
	      << " >> fx   : " << L2N_fx  << "\n"
	      << " >> fy   : " << L2N_fy  << "\n"
	      << " >> fx2  : " << L2N_fx2 << "\n"
	      << " >> fy2  : " << L2N_fy2 << "\n"
	      << " >> flap : " << L2N_lap << "\n"
		  ;

	write << "\nThe L inf Norm\n"
	      << " >> fx   : " << LInf_fx  << "\n"
	      << " >> fy   : " << LInf_fy  << "\n"
	      << " >> fx2  : " << LInf_fx2 << "\n"
	      << " >> fy2  : " << LInf_fy2 << "\n"
	      << " >> flap : " << LInf_lap << "\n"
		  ;

	write.close();

	std::cout << "[>] Done current step\n";

	return;
}