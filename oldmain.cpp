#ifndef INCLUDED_GLOBAL
#include "global.hpp"
#include <vector>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>   //tambahan ical
#include <cmath>     //tambahan ical
#include <chrono>    //tambahan ical
#endif

#ifndef INCLUDED_UTILS
#include "Utils.hpp"
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

// ===========================================================================
int main(int argc, char const *argv[])
{

#pragma region instances
	// TODO: creating instance of classes
	geometry geom_step;
	penalization penalization_step;
	poisson_solver advection_step;
	diffusion diffusion_step;
	remeshing remesh_step;
	save_data save_step;
	initialization initialization_step;
	LSMPSa lsmps;
	solverpoisson solver_poisson;
#pragma endregion

#pragma region internal_variables
	int nt_start = 0;
	double cum_time; // actual cumulative time
#pragma endregion

#pragma region variables
	// clock_t t1 = clock();
	Body body;		   // create obstacle object
	Particle particle; // create simulation particle object
	// t1 = clock() - t1;
	// printf("[%f]\n", (double)t1 / CLOCKS_PER_SEC);
#pragma endregion

	// ======================== DESCRIPTION ========================\
	------------------- start of computation --------------------
	printf("\n=============== Flow parameters ==================\n");
	printf("\tRe                = %12f \n", Pars::Re);
	printf("\tViscosity         = %12f \n", Pars::vis);
	printf("\tcore size         = %12f \n", Pars::sigma);
	printf("\tdt                = %12f \n", Pars::dt);
	printf("\tCourant number    = %12f \n", Pars::Courant);
	printf("\tTurbulent scaling = %12f \n", std::ceil(1.0e0 / Pars::Courant));
	printf("\tdelta time calculate   = %12f \n", 0.25*Pars::sigma*Pars::sigma/Pars::vis);
	printf("+------------------------------------------------+\n");

	// //=======================================================================================\
	// 	========================  IN CASE OF RUNNING VEM/EXTRA RESULTS ========================\
	// 	opt_extra_data = 1 or 2: running VEM code and without/with saving extra data\
	// 	opt_extra_data = 3: running Extra data calculation (by output)
	// if (Pars::opt_extra_data != 3) // running VEM
	// {

#pragma region geometry
	if(Pars::sim == 2 )
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
		// TODO: Geometry translational velocity
		geom_step.u_var(body); // generate the various traversal velocity	
	} else
	{
		geom_step.ball_sphere(body);
	}
#pragma endregion
//======================

#pragma region continue_simulation
	// =========================  New/"Continue the stopped" running  ========================\
		Only apply for moving body with constant Velocity\
		If deformable/accelerate moving body, we need save the vertex like particle's infomation
	if (Pars::opt_cont == 1)
	{
			
		// Open the file which have full data (before/at stop) , Current position of fluid-particles
		save_step.particle_data_reading(Pars::it_stop_full, particle.num, particle.x, particle.y, particle.s, particle.gz, particle.u, particle.v);//, Isumx, Isumy);
		
		// Current position of body-points
		for (int it = 0; it < Pars::it_stop_full; it++) // it_stop_full is changeable anyway, which depends on the last number of iteration when it is unfortunately stopped
		{
			for (int i = 0; i < body.x.size(); i++)
			{
				body.x[i] += body.uT[it] * Pars::dt;
				body.y[i] += body.vT[it] * Pars::dt;
			}
		}
		nt_start = Pars::it_stop_full;
	}
#pragma endregion

#pragma region restart_simulation
	else
	{
		if(Pars::sim == 2)
		{
			// TODO: Initial Particles Generation
			// initialization_step.init_particle(body, particle);
			// initialization_step.init_domain(particleBase);
			initialization_step.init_domain(particle);
			nt_start = 0; // because c-based indexing is started from '0'
		}else
		{
			initialization_step.init_domain_3d(particle);
			nt_start = 0;
		}
	}
#pragma endregion

	//=========================================================
	//TODO: SHOWS VIBRATION PARAMETERS
	if (Pars::vib == 1){
		
		Pars::mass = Pars::m_star * 0.5 * Pars::dens * pow(Pars::Df,2) - Pars::m_d;
		Pars::SpringConst = Pars::k_star * 0.5 * Pars::dens * pow(Pars::Uf,2);
		Pars::DamperConst = Pars::c_star * 0.5 * Pars::dens * Pars::Df * Pars:: Uf; 

		cout<<"Vibration Parameter"<<endl;
		printf("\tMass                = %12f \n", Pars::mass);
		printf("\tSpring Constant     = %12f \n", Pars::SpringConst);
		printf("\tDamper Constant     = %12f \n", Pars::DamperConst);
		printf("\tBody Area           = %12f \n", Pars::area);
		cout<<"---------------------------------"<<endl;

	} else if (Pars::vib == 2){
		
		Pars::inertia = Pars::i_star * 0.5 * Pars::dens * pow(Pars::lx,4); 
		Pars::SpringConst = pow(((Pars::uin / (Pars::U_star*Pars::lx)) * 2 * Pars::pi),2) * Pars::inertia;
		Pars::DamperConst = Pars::chi_star * 2 *sqrt(Pars::SpringConst*Pars::inertia); 
		Pars::tetha = Pars::alpha_a;
		
		cout<<"Vibration Parameter"<<endl;
		printf("\tInertia                = %12f \n", Pars::inertia);
		printf("\tSpring Constant     = %12f \n", Pars::SpringConst);
		printf("\tDamper Constant     = %12f \n", Pars::DamperConst);
		printf("\tAngle of Attack          = %12f \n", Pars::alpha_a);
		cout<<"---------------------------------"<<endl;
	}

	//CONTINUE TO THE SIMULATION?
	bool cek = false;
	string ans;
	if (cek == false){
		cout << "Continue? (yes/no)" << endl;
		cin>>ans;
		if (ans == "no"){
			cek = false;
		}else{
			cek = true; 
		}
	}

if (cek == true){
	// ===================================== LOOPS ===========================================
	cum_time = 0.0e0; //its cummulative time. atau waktu simulasi. 
	std::vector<double> _cumulativeTime;
	//clock_t delta_t; // actual time difference
	Pars::batas = -Pars::xdom + 5*Pars::sigma;
	
	for (size_t it = nt_start; it < Pars::nt; it++)
	{
		//delta_t = clock();
		auto t_start = std::chrono::system_clock::now();
		printf("+--------------- iter no. %d -------------------+\n", (int)it);
		// ? use displayBar

		//Karena partikel dipindah, maka butuh ngemesh ulang.
		// TODO: Remeshing
		if (it % (Pars::nrmsh) == 0)
		{
			remesh_step.get_remeshing(particle, body, it * Pars::dt, it);
		}
		cout << "succesfully remeshed" <<endl;
		//BOUNDARY 
		//DO NOT CHANGE BRUH
		if(it == 0)
		{
			for (size_t ii = 0; ii < particle.x.size();ii++){
				if (particle.x[ii] <= -Pars::xdom - 10 * Pars::sigma || particle.x[ii]  >= ((Pars::lxdom-Pars::xdom) + 6 * Pars::sigma)
					 || particle.y[ii]  <= (-Pars::lydom/2) - 10 * Pars::sigma || particle.y[ii]  >= ((Pars::lydom/2) + 6 * Pars::sigma ))
            	{
                	particle.isboundary.push_back(1);
					particle.boundaryval.push_back(0.0);
            	}
            	else
            	{
               		particle.isboundary.push_back(0);
					particle.boundaryval.push_back(0.0);
            	}
			}
		}

		if (it == 0){
			particle.vorticity.resize(particle.num, 0.0e0);
		}else{
			for (size_t i = 0; i < particle.num; i++)
			{
				particle.vorticity[i] = particle.gz[i] / pow(particle.s[i],2);
			}
		}

		if (it == 0){
			for (int j = 0; j < particle.x.size(); j++){
				particle.label.push_back(j);
			} 
		}

		// TODO: Poisson: solving Rotational Velocity, Stretching
		advection_step.poisson(particle, it);	
		//solver_poisson.set_solver(particle.u, particle.v, particle.x, particle.y, particle.s, particle.vorticity, particle.neighbor, particle.isboundary, particle.boundaryval, it);
		//lsmps.set_LSMPS_Laplace(particle.u, particle.x, particle.y, particle.s, particle.gz, particle.neighbor, particle.isboundary);


		// TODO: Correct vorticity and velocity using Brinkman penalization technique
		penalization_step.get_penalization(particle, body, it);

		//save_step.save_state(particle, "after_penalization_iterative");
		double test1,test2; 
		// TODO: Helmholtz decomposition
		// u = u_rotational + u_irrotational
		// where: u_rotational = u from Biot_Savart(vortex); u_irrotational = u_infinite
		// for iterative, calculate velocity again, then helmholtz decomposition
		// for classical, just do helmholtz decomposition
		if(Pars::iterative == 2){
			advection_step.poisson(particle, it);
			for (size_t i = 0; i < particle.num; i++)
			{
				particle.u[i] += Pars::uin;
				particle.v[i] += Pars::vin;
			}
			//d_base_save_data.force_pen(it, _p.num, lambda, kai, _p.u, _p.v, uS, vS, _p.s, xpus, _p.x, _p.y);
		} else {
			
			test1 = omp_get_wtime();
			for (size_t i = 0; i < particle.num; i++){
				particle.u[i] += Pars::uin;
				particle.v[i] += Pars::vin;
			}
			test2 = omp_get_wtime();
			printf("Velocity obtained (s) : %f s\n", test2-test1);
		}
		test1 = omp_get_wtime();
		save_step.Force2(it, 1,2,3,4, 1,2, particle);
		test2 = omp_get_wtime();
		printf("Force Calculation (Noca Method in s) : %f s\n", test2-test1);
		//if (it != 0) save_step.output(it, particle, body, cum_time);

		// TODO: Convection Sub-step
		std::vector<double> _dfdtConv(particle.num);
	
		advection_step.advection(particle, _dfdtConv); // ! later: do 2nd order scheme
		//advection_step.advection_rk2(particle, _dfdtConv);
		//selesai step 1
		
		//mulai step 2
		// TODO: Diffusion Sub-step
		std::vector<double> _dfdtDiff(particle.num);
		diffusion_step.main_diffusion(particle, _dfdtDiff); // ! later: do 2nd order scheme
		
		// time integration (1st order; Diffusion)
		for (size_t i = 0; i < particle.num; i++)
		{
			particle.gz[i] += Pars::dt * (_dfdtDiff[i]); 
		}
		
		
		//TODO: Move body (deformation, translation, rotation)
		geom_step.moving_body(it, body);
		//selesai step 2 dan ulangi
		
		// TODO: Saving data
		//delta_t = clock() - delta_t;							   // actual time [s]
		std::chrono::duration<double> elapsed_time_ms = (std::chrono::system_clock::now() - t_start);
		//_cumulativeTime.push_back(static_cast<double>(delta_t) / static_cast<double>(CLOCKS_PER_SEC));
		_cumulativeTime.push_back(elapsed_time_ms.count());
		cum_time = accumulate(_cumulativeTime.begin(), _cumulativeTime.end(), 0.0e0);
		// cum_time = cum_time + (double)delta_t / CLOCKS_PER_SEC; // cumulative time
		/*
		int _numberOfActiceParticle = 0;
		for (size_t i = 0; i < particle.num; i++)
		{
			if (particle.isActive[i])
			{
				_numberOfActiceParticle += 1;
			}
		}
		
		std::ofstream outs;
		if (it == 0)
		{
			outs.open("output/particle_number_time.dat");
			outs << it << "," << it * Pars::dt << "," << cum_time << "," << _numberOfActiceParticle << "\n";
			outs.close();
		}
		else if (it >= 1)
		{
			outs.open("output/particle_number_time.dat", std::ofstream::out | std::ofstream::app);
			outs << it << "," << it * Pars::dt << "," << cum_time << "," << _numberOfActiceParticle << "\n";
			outs.close();
		}
		std::cout << "\nnumber of particles: " << _numberOfActiceParticle << std::endl;*/
		std::cout << "\ncomputation time: " << cum_time << std::endl;
		//if (it % Pars::nt_sf == 0)
		//{
		save_step.output(it, particle, body, cum_time);

		//if (it == 1)
		//	break;
	}
}

	return 0;
}

// NOTEs:
// to avoid memory leak, please do an initalization for all struct members (e.g, set to be 0)

// std::ofstream outs;
// outs.open("output/data_temp.dat");
// for (size_t i = 0; i < particle.num; i++)
// 	outs << particle.x[i] << " " << particle.y[i] << " " << particle.gz[i] << "\n";
// outs.close();

// Body *b = new Body;
// std::vector<double> x = b->get_x();
// x.resize(1000000);
// b->set_x(x);
// cout << b->_x.size() << " " << x.size() << endl;

// Body *b = new Body;
// std::vector<double> &x = b->_x;
// x.resize(1000000);
// // b->set_x(x);
// cout << b->_x.size() << endl;

// delete b;