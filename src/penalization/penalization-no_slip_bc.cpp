#include "penalization.hpp"

// Simple penalization
void penalization::no_slip(Particle &_p, const Body &b)
{
	// Declaring internal variable
	// Derivative velocity
	std::vector<double> du(_p.num, 0.0e0);
	std::vector<double> dv(_p.num, 0.0e0);
	// Penalized velocity
	std::vector<double> u_pen(_p.num, 0.0e0);
	std::vector<double> v_pen(_p.num, 0.0e0);
	// Solid body velocity
	std::vector<double> uS(_p.num, 0.0e0);
	std::vector<double> vS(_p.num, 0.0e0);
	
	double uR, uT, uDEF, vR, vT, vDEF;

	// Calculate the solid velocity
	for (size_t i = 0; i < _p.num; i++)
	{
		// Rotation still not modeled
		uR = 0.0e0;      // Body node x-velocity influenced by body rotation
		vR = 0.0e0;      // Body node y-velocity influenced by body rotation

		// Tranlation take from the body variable
		uT = b.uT[0];    // Body node x-velocity influenced by body translation
		vT = b.vT[0];    // Body node y-velocity influenced by body translation
		
		// Still not consider the U_deformation yet!
		uDEF = 0.0e0;    // Body node x-velocity influenced by body deformation
		vDEF = 0.0e0;    // Body node y-velocity influenced by body deformation
		
		// Node resultant velocity
		uS[i] = uR + uT + uDEF;
		vS[i] = vR + vT + vDEF;
	}

	/*
	Penalization calculation procedure
	1. Calculate velocity penalization
	2. Update the vorticity
	3. Update the velocity and vorticity into the origin variable (Outside this function)
	Note:
	> There are implicit, semi-implicit, and explicit schemes
	> Velocity penalization calculation 1 for: Implicit and semi implicit
	> Velocity penalization calculation 2 for: Explicit
	*/

    // =================================================
	// =================== Process 1 ===================
	// =================================================
	// Velocity penalization calculation 1
	if (Pars::opt_pen == 1 || Pars::opt_pen == 2)
	{
		for (int i = 0; i < _p.num; i++)
		{
			// Consider all velocity
			u_pen[i] = (_p.u[i] + Pars::lambda*Pars::dt* _p.chi[i] * uS[i]) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));
			v_pen[i] = (_p.v[i] + Pars::lambda*Pars::dt* _p.chi[i] * vS[i]) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));

			// // Only consider the rotational velocity
			// u_pen[i] = (_p.u[i] + Pars::lambda * Pars::dt * _p.chi[i] * (uS[i] - Pars::u_inf)) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));
			// v_pen[i] = (_p.v[i] + Pars::lambda * Pars::dt * _p.chi[i] * (vS[i] - Pars::v_inf)) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));
		}
	}
	// Velocity penalization calculation 2
	else if (Pars::opt_pen == 3)
	{
		for (int i = 0; i < _p.num; i++)
		{
			// Consider all velocity
			u_pen[i] = (1 - _p.chi[i]) * _p.u[i] + _p.chi[i] * (uS[i]);
			v_pen[i] = (1 - _p.chi[i]) * _p.v[i] + _p.chi[i] * (vS[i]);
			
			// // Only consider the rotational velocity
			// u_pen[i] = (1 - _p.chi[i]) * _p.u[i] + _p.chi[i] * (uS[i] - Pars::u_inf);
			// v_pen[i] = (1 - _p.chi[i]) * _p.v[i] + _p.chi[i] * (vS[i] - Pars::v_inf);
		}
	}
	
	// Calculate the velocity difference
	for (size_t i = 0; i < _p.num; i++){
		du[i] = u_pen[i] - _p.u[i];
		dv[i] = v_pen[i] - _p.v[i];
	}
	
	// Update particle velocity
	for (int i = 0; i < _p.num; i++){
		_p.u[i] = u_pen[i];
		_p.v[i] = v_pen[i];
	}

	// =================================================
	// =================== Process 2 ===================
	// =================================================
	// Updating vorticity from calculating
	
	// Re-evaluate the neigbor (Note that in this region the resolution is single)
	_p.neighbor = this->d_neighbor.link_list(_p.num, _p.s, _p.x, _p.y, Pars::r_scale);
	
	// Performing LSMPS calculation
	LSMPSa lsmpsa_du;
	LSMPSa lsmpsa_dv;
	lsmpsa_du.set_LSMPS(_p.x, _p.y, _p.s, du, _p.neighbor);
	lsmpsa_dv.set_LSMPS(_p.x, _p.y, _p.s, dv, _p.neighbor);
	std::vector<double> _ddvdx = lsmpsa_dv.get_ddx();
	std::vector<double> _ddudy = lsmpsa_du.get_ddy();
	
	// Evaluating added penalization vorticity field [dw = grad (cross) du]
	double _dgpz, _dvor;
	for (int i = 0; i < _p.num; i++)
	{
		// Calculate the vorticity and vortex strength
		_dvor = _ddvdx[i] - _ddudy[i];         // w = del x u
		_dgpz = _dvor * std::pow(_p.s[i], 2);  // γ = w * Area

		// Update the vorticity and vortex strength
		_p.vorticity[i] += _dvor;
		_p.gz[i] += _dgpz;
	}

}

// // Iterative penalization
// void penalization::no_slip_iterative(Particle &_p, Particle &p, const Body &b, int it){
// 	int r_cut = Pars::r_scale; // ! neighbor_scale[DONE]

// 	std::vector<double> u0(_p.num, 0.0e0); //bikin vector, terus isinya langsung 0
// 	std::vector<double> v0(_p.num, 0.0e0); //bikin vector, terus isinya langsung 0
// 	std::vector<double> du(_p.num, 0.0e0); //bikin vector, terus isinya langsung 0
// 	std::vector<double> dv(_p.num, 0.0e0); //bikin vector, terus isinya langsung 0
// 	std::vector<double> uS(_p.num, 0.0e0); //Ini bagian solid !! 
// 	std::vector<double> vS(_p.num, 0.0e0); //Ini bagian solid !! 
// 	std::vector<double> dxi(_p.num, 0.0e0);
// 	std::vector<double> u_xi(_p.num, 0.0e0); 
// 	std::vector<double> v_xi(_p.num, 0.0e0); 
// 	std::vector<double> xi(_p.num, 0.0e0); 
	
// 	double uR, uT, uDEF, vR, vT, vDEF;
// 	//us masih case benda diam
// 	for (size_t i = 0; i < _p.num; i++)
// 	{
// 		uR = 0.0e0; 
// 		uT = b.uT[0]; 
// 		uDEF = 0.0e0;
// 		vR = 0.0e0; 
// 		vT = b.vT[0]; 
// 		vDEF = 0.0e0;
// 		uS[i] = uR + uT + uDEF;
// 		vS[i] = vR + vT + vDEF;
// 	}

// 	double Fx, Fy;
// 	Fx = 0.0; 
// 	Fy = 0.0;

// 	for (size_t i = 0; i < _p.num; i++){
// 		u0[i] = uS[i] - (_p.u[i]);
// 		v0[i] = vS[i] - (_p.v[i]);
// 	}

// 	_p.neighbor = d_neighbor.link_list(_p.num, _p.s, _p.x, _p.y, r_cut); 
	
// 	int nom = 0;
// 	for (int k = 1; k <= 20; k++){
// 		nom++;

// 		//calculate increment (alpha == 2)
// 		for(int i = 0; i < _p.num; i++){
// 			//du[i] = 2 * kai[i] * (u0[i] - u_xi[i]);
// 			//dv[i] = 2 * kai[i] * (v0[i] - v_xi[i]);
// 			du[i] = (u_xi[i] + lambda * Pars::dt * _p.chi[i] * (u0[i] - Pars::u_inf)) / (1.0e0 + (lambda * Pars::dt * _p.chi[i]));
// 			du[i] -= u_xi[i];
// 			dv[i] = (v_xi[i] + lambda * Pars::dt * _p.chi[i] * (v0[i] - Pars::v_inf)) / (1.0e0 + (lambda * Pars::dt * _p.chi[i]));
// 			dv[i] -= v_xi[i];
// 		}

// 		// ! GRADIENT HERE
// 		LSMPSa lsmpsa_du;
// 		LSMPSa lsmpsa_dv;
// 		lsmpsa_du.set_LSMPS(_p.x, _p.y, _p.s, du, _p.neighbor);
// 		lsmpsa_dv.set_LSMPS(_p.x, _p.y, _p.s, dv, _p.neighbor);
// 		std::vector<double> _ddvdx = lsmpsa_dv.get_ddx(); // [d(dv)/dx]
// 		std::vector<double> _ddudy = lsmpsa_du.get_ddy(); // [d(du)/dy]

// 		for (int i = 0; i < _p.num; i++)
// 		{
// 			dxi[i] = _ddvdx[i] - _ddudy[i];
// 			dxi[i] = 2 * (dxi[i]) * std::pow(_p.s[i], 2); //ini update strength
// 			xi[i] = xi[i] + dxi[i]; //ini jadi strength
// 		}

// 		//calculate new velocity using poisson
// 		int np = _p.num;
// 		double *x = new double[_p.num];
// 		double *y = new double[_p.num];
// 		double *s = new double[_p.num];
// 		double *gz = new double[_p.num];
// 		double *u = new double[_p.num];
// 		double *v = new double[_p.num];

// 		for (size_t i = 0; i < _p.num; i++)
// 		{	
// 			x[i] = _p.x[i];
// 			y[i] = _p.y[i];
// 			s[i] = _p.s[i];
// 			gz[i] = xi[i] * std::pow(_p.s[i] / Pars::sigma, 2); 
// 			u[i] = 0.0e0;
// 			v[i] = 0.0e0;
// 		}
// 		int _n0 = 1;
// 		int _iCutoff = Pars::icutoff;
// 		int _nS = Pars::n_s;
// 		int _nInter = 1;
// 		int _ndp = Pars::ndp;
// 		//printf(">> Iteration %d\n", k);
// 		FortranUtils::biotsavart_fmm_2d_(&_n0, &np, &np, x, y, s, u, v,
//           							 &_n0, &np, &np, x, y, s, gz,
//           							 &_iCutoff, &_nS, &_nInter, &_ndp);
		
// 		for (int i = 0; i < _p.num; i++)
// 		{
// 			u_xi[i] = u[i];
// 			v_xi[i] = v[i];
// 		}
// 		//FMM.biotsavart_fmm_2d( _particle, _particleDense,  _iCutoff, _nS,  _nInter,  _ndp);

// 		delete[] x;
// 		delete[] y;
// 		delete[] s;
// 		delete[] gz;
// 		delete[] u;
// 		delete[] v;

// 		//calculate added penalization force //check lagi
// 		double dfx = 0.0, dfy = 0.0;
// 		for (size_t i = 0; i < _p.num; i++){
// 			double _area = std::pow(_p.s[i] / Pars::sigma, 2);
//         	dfx += _p.y[i] * dxi[i] * _area;
//        		dfy += -_p.x[i] * dxi[i] * _area;
// 		}

// 		dfx = (-Pars::RHO / Pars::dt) * dfx;
// 		dfy = (-Pars::RHO / Pars::dt) * dfy;
// 		//printf("%f, %f\n",dfx, dfy);
// 		Fx = Fx + dfx;
// 		Fy = Fy + dfy;
	
// 		//break condition (epsilon = 10^-1)
// 		if (std::sqrt(std::pow(dfx,2) + std::pow(dfy,2)) /  std::sqrt(std::pow(Fx,2) + std::pow(Fy,2)) < 0.01){
// 			break;
// 		}
// 	}

// 	for (int i = 0; i < _p.num; i++){
// 		_p.gz[i] += dxi[i]; // add the iterative vorticity into strength by multiplying it with size of particle
// 	}

// 	double xpus[2];
// 	// COMMENTED <========= UNCOMMENT IF NECCESARRY
// 	// d_base_save_data.force_pen(it, _p.num, lambda, kai, _p.u, _p.v, uS, vS, _p.s, xpus, _p.x, _p.y);

// 	printf(">>Finished iteration in %d\n", nom);
// 	double cd = 2 * Fx / (Pars::RHO * Pars::u_inf * Pars::u_inf * 1); //area masih 1
// 	double cl = 2 * Fy / (Pars::RHO * Pars::u_inf * Pars::u_inf * 1);

// 	//save data
// 	std::ofstream ofs;
//     if(it == 0){
//         printf("v_inf force from iterative penalization....\n");
//         ofs.open("force-iterative.csv");
//         ofs <<"Iter" << "," << "Cd" << "," << "Cl" << "," <<"Fx" 
//             << "," << "Fy\n";
//         ofs << it << "," << cd << "," << cl << "," << Fx
//             << "," << Fy << "\n";
//         ofs.close();
//     }else{
//         printf("v_inf force from iterative penalization....\n");
//         ofs.open("force-iterative.csv", std::ofstream::out | std::ofstream::app);
// 	    ofs << it * Pars::dt << "," << cd << "," << cl << "," << Fx
//             << "," << Fy << "\n";
// 	    ofs.close();
//     }

// 	//check again;
// 	for(int s= 0; s<_p.num;s++){
// 		_p.gz[s] = _p.gz[s] + xi[s];
// 	}
// }