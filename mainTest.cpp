#include "Utils.hpp"

#ifndef INCLUDED_GLOBAL
#include "global.hpp"
#include <vector>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>   // Ical added
#include <cmath>     // Ical added
#include <chrono>    // Ical added
#endif

#ifndef INCLUDED_GEOMETRY
#include "src/geometry/geometry.hpp"
#endif

#ifndef INCLUDED_POISSON_SOLVER
#include "src/advection/poisson.hpp"
#endif

#ifndef INCLUDED_DC_OPERATOR
#include "src/diffusion/diffusion.hpp"
#endif

#ifndef INCLUDED_REMESHING
#include "src/remeshing/remeshing.hpp"
#endif

#ifndef INCLUDED_PENALIZATION
#include "src/penalization/penalization.hpp"
#endif

#ifndef INCLUDED_SAVE_DATA
#include "src/save_data/save_data.hpp"
// Data stored at each or some iteration:
	// The simulation parameter (number, domain size, etc.) 	>>> make a log
	// Run time counter (at each iteration) 		    		>>> Data list per count iteration
	// Particle data (Neighbor list, chi, active particle, etc.)>>> Data list per count iteration
	// Property contour data (V, w, P, T)						>>> Data list per count iteration
	// Additional data (drag, force, moment, lift)				>>> Data list per count iteration
#endif

#ifndef INCLUDED_LSMPS
#include "src/LSMPS/LSMPSa.hpp"
#endif

#ifndef INCLUDED_INITIALIZATION
#include "src/initialization/initialization.hpp"
#endif

#ifndef INCLUDED_FORTRAN_UTILS
#include "src/Fortran/Utils/FortranUtils.hpp"
#endif

#ifndef INCLUDED_DC_BASE
#include "src/DC_operator/base_dc.hpp"
#endif

#ifndef INCLUDED_SOLVER_POISSON
#include "src/solver_poisson/solver.hpp"
#endif

#ifndef INCLUDED_SAVE_DATA_BASE
#include "src/save_data/base_save_data.hpp"
#endif

// ===========================================================================
int main(int argc, char const *argv[])
{

#pragma region instances 	// TODO: creating instance of subroutine
	geometry geom_step;
	penalization penalization_step;
	poisson_solver advection_step;
	diffusion diffusion_step;
	remeshing remesh_step;
	save_data save_step;
	initialization initialization_step;
	LSMPSa lsmps;
	solverpoisson solver_poisson;
	base_save_data d_base_save_data;
#pragma endregion

#pragma region variables	// TODO: creating instance of object
	Body body;		   		// The solid object
	Particle particle; 		// The simulation particle object
#pragma endregion

#pragma region internal_variables	// This is used for continuing simulation
	int nt_start = 0;	// The time start
	double cum_time; 	// actual cumulative time
#pragma endregion

	// Initial Parameter Summary Data
	// Shoud be moved in different program
	printf("\n=============== Flow parameters ==================\n");
	printf("\tRe                = %12f \n", Pars::RE);
	printf("\tViscosity         = %12f \n", Pars::NU);
	printf("\tcore size         = %12f \n", Pars::sigma);
	printf("\tdt                = %12f \n", Pars::dt);
	printf("\tCourant number    = %12f \n", Pars::Courant);
	printf("\tCourant2 (<0.5)   = %12f \n", Pars::Diffusion);
	printf("\tTurbulent scaling = %12f \n", std::ceil(1.0e0 / Pars::Courant));
	printf("\tdelta time calculate   = %12f \n", 0.25*Pars::sigma*Pars::sigma/Pars::NU);
	printf("+------------------------------------------------+\n");

	// Saving the simulation parameter
	std::ofstream param;
	param.open("output/Parameter.dat");
	param << "Simulation Flow Parameter:\n";
	param << "> stream_velocity : " << Pars::u_inf  << " m/s"<< "\n"
	      << "> viscosity    : " << Pars::NU  << " m^2/s" << "\n"
		  << "> density      : " << Pars::RHO << " kg/m^3" << "\n"
		  << "> plate_length : " << Pars::lx   << " m" << "\n"
		  << "> plate_thick  : " << Pars::lx*Pars::H_star   << " m" << "\n"
		  << "> Re           : " << Pars::RE   << " [-]" << "\n"
		  << "> diameter     : " << Pars::sigma<< " m" << "\n"
		  << "> time_step    : " << Pars::dt   << " s" << "\n"
		  << "> xdom_length  : " << Pars::lxdom   << " m" << "\n"
		  << "> ydom_length  : " << Pars::lydom   << " m" << "\n"
		  << "> plate pos    : " << Pars::xdom   << " m" << "\n"
		  << "> sim_time     : " << Pars::simulation_time   << " s" << "\n";
		  
	// param.close();
	int particle_max_number = 0;


	// //===========================================================================\
	// 	=============  IN CASE OF RUNNING VEM/EXTRA RESULTS ========================\
	// 	opt_extra_data = 1 or 2: running VEM code and without/with saving extra data\
	// 	opt_extra_data = 3: running Extra data calculation (by output)

#pragma region geometry
	if(Pars::DIM == 2 )
	{
		// TODO: Geometry
		if (Pars::opt_body == 1)
			geom_step.cylinder_generator(body); // generates the cylinder
		else if (Pars::opt_body == 2)
			geom_step.naca_generator(body); // generates the NACA 4digit
		else if (Pars::opt_body == 3)
			geom_step.plate_generator(body); // generates the plate
		else if (Pars::opt_body == 4)
			geom_step.square_generator(body);
		else if (Pars::opt_body == 5)
		 	geom_step.flat_plate_generator(body);
		// else if (Pars::opt_body == 6)
		//  	geom_step.flat_wedge_plate_generator(body);
		// TODO: Geometry translational velocity
		geom_step.u_var(body);	// generate the various traversal velocity	
	} else
	{
		geom_step.ball_sphere(body);
	}
#pragma endregion

	// ============= Initialization
	// TODO: Initial Particles Generation
	if(Pars::DIM == 2)
	{
		// initialization_step.init_particle(body, particle);	// <<<========= Type 1
		// initialization_step.init_domain(particleBase);		// <<<========= Type 2
		initialization_step.init_domain_2d(particle);			// <<<========= Type 3
	}else
	{
		initialization_step.init_domain_3d(particle);
	}

	//=========================================================
	//TODO: SHOWS VIBRATION PARAMETERS
	if (Pars::vib == 1){
		
		Pars::mass = Pars::m_star * 0.5 * Pars::RHO * pow(Pars::Df,2) - Pars::m_d;
		Pars::SpringConst = Pars::k_star * 0.5 * Pars::RHO * pow(Pars::U_inf,2);
		Pars::DamperConst = Pars::c_star * 0.5 * Pars::RHO * Pars::Df * Pars:: U_inf; 

		cout<<"Vibration Parameter"<<endl;
		printf("\tMass                = %12f \n", Pars::mass);
		printf("\tSpring Constant     = %12f \n", Pars::SpringConst);
		printf("\tDamper Constant     = %12f \n", Pars::DamperConst);
		printf("\tBody Area           = %12f \n", Pars::area);
		cout<<"---------------------------------"<<endl;

	} else if (Pars::vib == 2){
		
		Pars::inertia = Pars::i_star * 0.5 * Pars::RHO * pow(Pars::lx,4); 
		Pars::SpringConst = pow(((Pars::u_inf / (Pars::U_star*Pars::lx)) * 2 * Pars::pi),2) * Pars::inertia;
		Pars::DamperConst = Pars::chi_star * 2 *sqrt(Pars::SpringConst*Pars::inertia); 
		Pars::tetha = Pars::alpha_a;
		
		cout<<"Vibration Parameter"<<endl;
		printf("\tInertia                = %12f \n", Pars::inertia);
		printf("\tSpring Constant     = %12f \n", Pars::SpringConst);
		printf("\tDamper Constant     = %12f \n", Pars::DamperConst);
		printf("\tAngle of Attack          = %12f \n", Pars::alpha_a);
		cout<<"---------------------------------"<<endl;
	}

	// //CONTINUE TO THE SIMULATION?
	// bool cek = false;
	// string ans;
	// if (cek == false){
	// 	cout << "Continue the simulation? (yes/no)" << endl;
	// 	cin >> ans;
	// 	if (ans == "no"){
	// 		cek = false;
	// 	}else{
	// 		cek = true; 
	// 	}
	// }

// if (cek == true){
// 	// ===================================== LOOPS ===========================================
// 	cum_time = 0.0e0; 						// The cummulative time through the simulation
// 	std::vector<double> _cumulativeTime;	// Store the cummulative time of each iteration
// 	//clock_t delta_t; // actual time difference
// 	// Pars::batas = -Pars::xdom + 5*Pars::sigma;

// 	for (size_t it = nt_start; it < Pars::nt; it++)
// 	{
// 		//delta_t = clock();
// 		auto t_start = std::chrono::system_clock::now();
// 		printf("+--------------- iter no. %d -------------------+\n", (int)it);

// 		// TODO: Remeshing
// 		if (it % (Pars::nrmsh) == 0) // remesh every nrmsh time step
// 		{
// 			remesh_step.get_remeshing(particle, body, it * Pars::dt, it);
// 		}
// 		cout << "succesfully remeshed" <<endl;

// 		// BOUNDARY (Is this really is important ??)
// 		// PLEASE DO NOT CHANGE
// 		if(it == 0)
// 		{
// 			for (size_t ii = 0; ii < particle.x.size();ii++){
// 				if (particle.x[ii] <= -Pars::xdom - 7 * Pars::sigma || particle.x[ii]  >= ((Pars::lxdom-Pars::xdom) + 6 * Pars::sigma)
// 					 || particle.y[ii]  <= (-Pars::lydom/2) - 7 * Pars::sigma || particle.y[ii]  >= ((Pars::lydom/2) + 6 * Pars::sigma ))
//             	{
//                 	particle.isboundary.push_back(1);
// 					particle.boundaryval.push_back(0.0);
//             	}
//             	else
//             	{
//                		particle.isboundary.push_back(0);
// 					particle.boundaryval.push_back(0.0);
//             	}
// 			}
// 		}

// 		if (it == 0){
// 			particle.vorticity.resize(particle.num, 0.0e0);
// 		}else{
// 			for (size_t i = 0; i < particle.num; i++)
// 			{
// 				particle.vorticity[i] = particle.gz[i] / pow(particle.s[i],2);	// Why not the particle area (?)
// 			}
// 		}

// 		if (it == 0){
// 			for (int j = 0; j < particle.x.size(); j++){
// 				particle.label.push_back(j);
// 			}
// 		}

// 		// ---- COMMENT ALL THE CODE BELOW FOR INITIAL CHECK ----//

// 		// TODO: Poisson: solving Rotational Velocity, Stretching
// 		advection_step.poisson(particle, it);
// 		//solver_poisson.set_solver(particle.u, particle.v, particle.x, particle.y, particle.s, particle.vorticity, particle.neighbor, particle.isboundary, particle.boundaryval, it);
// 		//lsmps.set_LSMPS_Laplace(particle.u, particle.x, particle.y, particle.s, particle.gz, particle.neighbor, particle.isboundary);

// 		// TODO: Correct vorticity and velocity using Brinkman penalization technique
// 		penalization_step.get_penalization(particle, body, it);

// 		//save_step.save_state(particle, "after_penalization_iterative");
// 		double test1,test2; 
// 		// TODO: Helmholtz decomposition
// 		// u = u_rotational + u_irrotational
// 		// where: u_rotational = u from Biot_Savart(vortex); u_irrotational = u_infinite
// 		// for iterative, calculate velocity again, then helmholtz decomposition
// 		// for classical, just do helmholtz decomposition
// 		if(Pars::iterative == 2){
// 			advection_step.poisson(particle, it);
// 			for (size_t i = 0; i < particle.num; i++)
// 			{
// 				particle.u[i] += Pars::uin;
// 				particle.v[i] += Pars::vin;
// 			}
// 			//d_base_save_data.force_pen(it, _p.num, lambda, kai, _p.u, _p.v, uS, vS, _p.s, xpus, _p.x, _p.y);
// 		} else {
			
// 			test1 = omp_get_wtime();
// 			for (size_t i = 0; i < particle.num; i++){
// 				particle.u[i] += Pars::uin;
// 				particle.v[i] += Pars::vin;
// 			}
// 			test2 = omp_get_wtime();
// 			printf("Velocity obtained (s) : %f s\n", test2-test1);
// 		}
		
// 		// Force calculation and saving  <<<<< Added by Angga 18 May 2022
// 		//test1 = omp_get_wtime();
// 		if (Pars::force_type == 1){
// 			double xpus[2] = {0, 0};
// 			d_base_save_data.force_pen(it, particle, xpus);
// 		}
// 		else if (Pars::force_type == 2){
// 			save_step.Force2(it, 1,2,3,4, 1,2, particle);
// 		}
		
// 		// TODO: Convection Sub-step
// 		std::vector<double> _dfdtConv(particle.num);
	
// 		advection_step.advection(particle, _dfdtConv); // ! later: do 2nd order scheme
		
// 		//mulai step 2
// 		// TODO: Diffusion Sub-step
// 		std::vector<double> _dfdtDiff(particle.num);
// 		diffusion_step.main_diffusion(particle, _dfdtDiff); // ! later: do 2nd order scheme
		
// 		// time integration (1st order; Diffusion)
// 		for (size_t i = 0; i < particle.num; i++)
// 		{
// 			particle.gz[i] += Pars::dt * (_dfdtDiff[i]); 
// 		}
		
// 		//TODO: Move body (deformation, translation, rotation)
// 		// geom_step.moving_body(it, body);						// Unccoment for moving body <======
// 		//selesai step 2 dan ulangi
		
// 		// TODO: Saving data
// 		//delta_t = clock() - delta_t;							   // actual time [s]
// 		std::chrono::duration<double> elapsed_time_ms = (std::chrono::system_clock::now() - t_start);
// 		//_cumulativeTime.push_back(static_cast<double>(delta_t) / static_cast<double>(CLOCKS_PER_SEC));
// 		_cumulativeTime.push_back(elapsed_time_ms.count());
// 		cum_time = std::accumulate(_cumulativeTime.begin(), _cumulativeTime.end(), 0.0e0);
// 		// cum_time = cum_time + (double)delta_t / CLOCKS_PER_SEC; // cumulative time
// 		/*
// 		int _numberOfActiceParticle = 0;
// 		for (size_t i = 0; i < particle.num; i++)
// 		{
// 			if (particle.isActive[i])
// 			{
// 				_numberOfActiceParticle += 1;
// 			}
// 		}
		
// 		std::ofstream outs;
// 		if (it == 0)
// 		{
// 			outs.open("output/particle_number_time.dat");
// 			outs << it << "," << it * Pars::dt << "," << cum_time << "," << _numberOfActiceParticle << "\n";
// 			outs.close();
// 		}
// 		else if (it >= 1)
// 		{
// 			outs.open("output/particle_number_time.dat", std::ofstream::out | std::ofstream::app);
// 			outs << it << "," << it * Pars::dt << "," << cum_time << "," << _numberOfActiceParticle << "\n";
// 			outs.close();
// 		}
// 		std::cout << "\nnumber of particles: " << _numberOfActiceParticle << std::endl;*/
// 		std::cout << "\ncomputation time: " << cum_time << std::endl;
// 		//if (it % Pars::nt_sf == 0)
// 		//{
				
// 		// save_step.output(it, particle, body, cum_time);				// Unccoment for moving body <======

// 		// Save the particle data for each step after remeshing
// 		int time_temp = std::ceil(Pars::nt/Pars::nt_data);
// 		if (it % time_temp == 0){
// 			int addDigit = std::floor(std::log10(Pars::nt)) - std::floor(std::log10(it+1));
// 			std::string DataName;
// 			for (int digit = 0; digit < addDigit; digit++)
// 				DataName.append("0");
// 			DataName.append(to_string(it));
// 			save_step.save_state(particle, DataName, 2);		// New saving particle data
// 		}

// 		// if (it == 2) break;
		
// 		particle_max_number = particle.num > particle_max_number ? particle.num : particle_max_number;
// 	}
	
// 	param << "> particle_number : " << particle_max_number << "\n";
// 	param.close();
	
// 	// Saving the simulation time for each iteration
// 	std::ofstream time_note;
// 	double timer_data;
// 	time_note.open("output/Simulation_Time.dat");
// 	time_note << "iteration,sim_time,comp_time\n";
// 	for (int num = 0; num < _cumulativeTime.size(); num++){
// 		timer_data = num * Pars::dt;
// 		time_note << "" << num << "," << timer_data << "," << _cumulativeTime[num] << "\n";
// 	}
// 	time_note << "\nTotal Simulation time = " << cum_time << "\n";
// 	time_note.close();
// }
// 	std::cout << "There are " << argc << " arguments:" << std::endl;
// 	for (size_t i = 0; i < argc; i ++){
// 		std::cout << argv[i] << std::endl;
// 	}

	return 0;
}

// NOTEs:
// to avoid memory leak, please do an initalization for all struct members (e.g, set to be 0)

// The simulation outline:
// 1. Set the parameter
//    > The domain size
//    > Particle size
//    > Simulation time
//    > All simulation parameter
// 2. Initialization (
// 3. Run all code (Follow the initial flow code)



// c:\program files\gnuwin32\bin\make.exe
// Test!
// 1. Check the initialization subroutine
// 2. Check the remeshing subroutine
// >>> Need to make new output file
// 3. Start the full simulation (coarse >then> fine)
