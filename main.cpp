#pragma region included_package
	#ifndef INCLUDED_UTILS		// Already include global
	#include "Utils.hpp"
	#endif

	// Subroutine Packages
	// *******************
	#ifndef INCLUDED_GEOMETRY
	#include "src/geometry/geometry.hpp"
	#endif

	#ifndef INCLUDED_INITIALIZATION
	#include "src/initialization/initialization.hpp"
	#endif

	#ifndef INCLUDED_REMESHING
	#include "src/remeshing/remeshing.hpp"
	#endif

	#ifndef INCLUDED_DATA_SAVING
	#include "src/save_data/data_saving.hpp"
	#endif

	// Physics Calculation Method
	// **************************
	#ifndef INCLUDED_PENALIZATION
	#include "src/penalization/penalization.hpp"
	#endif

	#ifndef INCLUDED_ADVECTION
	#include "src/advection/advection.hpp"
	#endif

	#ifndef INCLUDED_DIFFUSION
	#include "src/diffusion/diffusion.hpp"
	#endif

	#ifndef INCLUDED_VELOCITY_POISSON
	#include "src/velocity_poisson/velocity_poisson.hpp"
	#endif

	#ifndef INCLUDED_PRESSURE_POISSON
	#include "src/pressure_poisson/pressure_poisson.hpp"
	#endif

	void saveResidual(Particle&, int);

#pragma endregion

// ************* Program Starts Here ************* //
int main(int argc, char const *argv[])
{
	#pragma region data_storage
		Body body;		    // The solid object
		Particle particle; 	// The simulation particle object
	#pragma endregion

	#pragma region subroutine_instances
		// Initialization tools
		geometry geom_step;                  // Geometry generation
		initialization initialization_step;  // Particle distribution
		remeshing remesh_step;               // Particle redistribution and neighbor search
		
		// Solver tools
		penalization penalization_step;      // Penalization
		velocity_poisson velocity_step;      // FMM:Biot Savart velocity solver
		advection advection_step;       	 // Advection
		diffusion diffusion_step;            // Diffusion
		// solverpoisson solver_poisson;        // Poisson solver
		pressure_poisson pressure_step;		 // Poisson solver of pressure
		
		// Simulation tools
		simUtil utilitis_step;               // Simulation utilities
		save_data save_step;                 // Save particle distribution
		base_save_data d_base_save_data;     // Save force each time interval
		std::ofstream _data;
	#pragma endregion

	// [LOG] Display Simulation Log
	save_step.summary_log();

	#pragma region internal_variable
		double curr_comp_time = 0.0e0;  	// Total computational time at each iteration
		double cum_comp_time = 0.0e0;   	// The total cumulative of simulation computational time
		int nt_start = Pars::cont_num;	    // The starting iteration step
	#pragma endregion

	#pragma region stability_variable
		double CourMax = Pars::Courant;	    // The global maximum courant number throughout simulation
		double DiffMax = Pars::Diffusion;	// The global maximum diffusion number throughout simulation
		double StabMax = 50*Pars::sigma*Pars::sigma/Pars::NU;	// The global maximum stability criteria throughout simulation
	#pragma endregion
	
	// Pre processing prompt displaying
	printf("\n#=================================================#");
	printf("\n+---------------- PRE PROCESSING -----------------+");
	printf("\n#=================================================#\n");

	// ========== Body Region ==========
	// *********************************
	geom_step.generateBody(body);			// Generate the body data
	save_step.save_state(body,"initial");	// Save the body data

	// ========== Particle Region ==========
	// *************************************
	// Initial particle distribution generation
	if (Pars::opt_sim_cont == 0){				// Start at initial time
		initialization_step.generate_initial(body, particle);
		nt_start = 0;
	}else if(Pars::opt_sim_cont == 1){			// Continue simulation (from state data)
		initialization_step.generate_continue(body, particle, nt_start);
		nt_start += 1;
	}
	
	// Generate cell list and perform neighbor search
	remesh_step.remeshing_init(particle);
	
	// Save the particle data
	if (Pars::opt_sim_cont == 0){				// Start at initial time
		save_step.save_state(particle,"initial", 0);
	}else if (Pars::opt_sim_cont == 1){			// Continue simulation (from state data)
		save_step.save_state(particle,"continue", 0);
	}

	/*
	// ========== Vibrational Parameter Display ==========	
	if (Pars::vib == 1){
		
		Pars::mass = Pars::m_star * 0.5 * Pars::RHO * pow(Pars::Df,2) - Pars::m_d;
		Pars::SpringConst = Pars::k_star * 0.5 * Pars::RHO * pow(Pars::U_inf,2);
		Pars::DamperConst = Pars::c_star * 0.5 * Pars::RHO * Pars::Df * Pars:: U_inf; 

		printf("\n+-------------- Vibration Parameter --------------+\n");
		printf("Mass                                    : %8.4f \n", Pars::mass);
		printf("Spring Constant                         : %8.4f \n", Pars::SpringConst);
		printf("Damper Constant                         : %8.4f \n", Pars::DamperConst);
		printf("Body Area                               : %8.4f \n", Pars::area);
		printf("+-------------------------------------------------+\n");

	} else if (Pars::vib == 2){
		
		Pars::inertia = Pars::i_star * 0.5 * Pars::RHO * pow(Pars::Df,4); 
		Pars::SpringConst = pow(((Pars::u_inf / (Pars::U_star*Pars::Df)) * 2 * Pars::pi),2) * Pars::inertia;
		Pars::DamperConst = Pars::chi_star * 2 *sqrt(Pars::SpringConst*Pars::inertia); 
		Pars::tetha = Pars::alpha_a;
		
		printf("\n+-------------- Vibration Parameter --------------+\n");
		printf("Inertia                                 : %8.4f \n", Pars::inertia);
		printf("Spring Constant                         : %8.4f \n", Pars::SpringConst);
		printf("Damper Constant                         : %8.4f \n", Pars::DamperConst);
		printf("Angle of Attack                         : %8.4f \n", Pars::alpha_a);
		printf("+-------------------------------------------------+\n");
	}
	*/
	
	// ========== Initialization Summary ==========
	// ********************************************
	printf("\n+------------ Initialization Summary -------------+\n");
	printf("Number of body node                     : %8d \n", body.num);
	printf("Number of particle node                 : %8d \n", particle.num);
	printf("Total iteration number                  : %8d \n", Pars::nt);
	printf("+-------------------------------------------------+\n");
	
	// [LOG] Save Summary Log
	save_step.save_summary_log(particle);

	// ========== Simulation Run Prompt ==========
	// *******************************************
	// Simulation command prompt
	std::cout << std::endl <<
		"<!> The initialization is completed!\n" << 
		"<!> Do you want to run the simulation? (yes/no)\n" <<
		"    Type here: ";
	// Prompt to continue run the simulation
	bool _run = false;std::string _cmd;std::cin >> _cmd;
	if (_cmd == "yes" | _cmd == "Yes")	// RUN
	{_run = true;}
	else								// DON'T RUN
	{_run = false;}

	// =====================================================
	// =============== SIMULATION ITERATION ================
	// =====================================================
	if (_run == true)
	{
		// Solving processing header prompt
		printf("\n#=================================================#");
		printf("\n+-------------------- SOLVING --------------------+");
		printf("\n#=================================================#");		
		bool resetFMMTree = false;
		
		// *************** Iteration Loop Starts Here *************** //
		for (int step = nt_start; step < Pars::nt; step++)
		{
			// *************** ITERATION INITIAL SETTING and DISPLAY *************** //
			// Calculate the number of iteration digit
			utilitis_step.startCounter(step);
			utilitis_step.printHeader(step);

			// Printing the stability criteria: courant (C) and diffusion (Phi) number
			if (Pars::flag_disp_stability){
				std::vector<double> max_stab;
				utilitis_step.stabilityEval(particle, max_stab);

				// Update the global stability criteria maximum value
				if (max_stab[0] > CourMax){CourMax = max_stab[0];}
				if (max_stab[1] > DiffMax){DiffMax = max_stab[1];}
				if (max_stab[2] > StabMax){StabMax = max_stab[2];}
			}

			// Solver computational time manager
			// clock_t compTime = clock();		// Using time_t
			auto t_start = std::chrono::system_clock::now();	// Using chrono

			// *************** THE SOLVER *************** //
			{
				/*
				Sequence of solving computation:
				1. Particle Redistribution  -> Structured particle distribution
				2. Velocity calculation     -> Structured particle distribution
				   -> Calculating rotational velocity
				   -> Calculating total velocity by helmholtz decomposition
				   -> Saving particle step [The best stage to save particle]
				3. Penalization             -> Structured particle distribution
				   -> Saving particle force [Must be save after penalization]
				   -> Saving particle step [Alternative stage]
				4. Perform Advection        -> Unstructured particle distribution 
				5. Perform Diffusion        -> Unstructured particle distribution
				   -> Saving particle step [Unstructured distribution]
				
				Try to rearrange the solver seq. 3->4->5->1->2
				*/
				
				// SOLVER SEQUENCE : 3 -> 4 -> 5 -> 1 -> 2

				// =============== SOLVER STEP [3] ===============
				// [3] Perform penalization using Brinkmann: Penalize the velocity in body domain
				int __step;
				if ((Pars::opt_sim_cont == 1) && (step == nt_start + 1)){__step = 0;}else{__step = step;}
				
				penalization_step.get_penalization(particle, body, __step);
				
				// Iterative penalization calculation
				if(Pars::opt_pen_iter == 2){
					// Calculate the rotational velocity
					velocity_step.get_velocity(particle, __step);
					// Update the velocity by helmholtz decomposition
					for (size_t i = 0; i < particle.num; i++)
					{
						particle.u[i] += Pars::u_inf;
						particle.v[i] += Pars::v_inf;
					}
					
					// Calculate the penalization
					penalization_step.get_penalization(particle, body, __step);
				}

				// ================ Barrier Mark =================
				// [!] TODO: Saving Data Force 
				// Force calculation and saving [Inserted after the Penalization calculation]
				std::cout << "\nSaving force data ...\n";
				if (Pars::opt_force_type == 1){
					std::cout << "<+> Penalization Force\n";
					double xpus[2] = {0, 0};
					d_base_save_data.force_pen(step, particle, xpus);
				}
				else if (Pars::opt_force_type == 2){
					save_step.Force2(step,1,2,3,4,1,2,particle);
				}
				
				// =============== SOLVER STEP [4] ===============
				// [4] Convection/Advection Sub-step: Perform the particle advection
				advection_step.advection_euler(particle);      // ! later: do 2nd order scheme
				
				// =============== SOLVER STEP [5] ===============
				// [5] Diffusion Sub-step: Calculate the vorticity diffusion & time integration
				diffusion_step.main_diffusion(particle); // ! later: do 2nd order scheme

				// =============== SOLVER STEP [1] ===============
				// [1] Particle redistribution: rearrange the particle distribution by interpolating vorticity
				if ((step % Pars::nrmsh) == 0) // redistribute particle every given iteration step
				{
					resetFMMTree = remesh_step.get_remeshing(particle, body, step);
				}

				// =============== SOLVER STEP [2] ===============
				// [2] Velocity calculation by poisson/biot savart: solving Rotational Velocity & Stretching
				// int __step;
				if (((Pars::opt_sim_cont == 1) && (step == nt_start + 1)) || resetFMMTree){__step = 0;}else{__step = step;}
				// std::cout << "The value of _step : " << __step << "\n";
				velocity_step.get_velocity(particle, __step);
				
				// Helmholtz decomposition
				// u = u_rotational + u_irrotational
				printf("<+> Adding the irrotational velocity term \n");
				printf("    [helmholtz backward decomposition]\n");
				for (size_t i = 0; i < particle.num; i++){
					particle.u[i] += Pars::u_inf;
					particle.v[i] += Pars::v_inf;
				}

			}	// End of the solver calculation

			// *************** POST PROCESSING *************** //
			// Save the particle data at given step interval
			if (step % Pars::step_inv == 0){
				// Pressure calculation
				if (Pars::flag_pressure_calc){
					pressure_step.get_pressure(particle);
				}
				
				// Saving particle data
				std::string DataName;
				utilitis_step.saveName(DataName, step);			// Writing data file name
				save_step.save_state(particle, DataName, 0);	// Saving particle data
				
				// Saving Residual
				if (Pars::flag_save_residual == true){
					saveResidual(particle,step);
				}
			}

			// *************** COMPUTATIONAL TIME *************** //
			// Calculate the accumulative computational time
			std::chrono::duration<double> elapsed_time_ms = (std::chrono::system_clock::now() - t_start);
			curr_comp_time = elapsed_time_ms.count();				// Current iteration computational time
			cum_comp_time += curr_comp_time;						// Accumulative computational time

			// Saving the simulation time for each iteration
			if (Pars::flag_save_sim_time){
				if (step == 0){
					_data.open("output/Simulation_Time.csv");
					_data << "iteration,sim_time,par_num,comp_time,curr_cum_comp_time\n";
					_data << step << "," << step * Pars::dt << "," << particle.num << "," << curr_comp_time << "," << cum_comp_time << "\n";
					_data.close();
				}else{
					_data.open("output/Simulation_Time.csv", std::ofstream::out | std::ofstream::app);
					_data << step << "," << step * Pars::dt << "," << particle.num << "," << curr_comp_time << "," << cum_comp_time << "\n";
					_data.close();
				}
			}

			// Displaying the simulation time for each iteration
			printf("\n**************** Iteration Summary ****************");
			printf("\n<!> Current iteration sim. time      : %9.3f s", step * Pars::dt);
			printf("\n<!> Current iteration comp. time     : %9.3f s", curr_comp_time);
			printf("\n<!> Cumulative computational time    : %9.3f s\n", cum_comp_time);

			// Prediction time to finish the simulation
			if (Pars::flag_disp_pred_time == true){
				utilitis_step.predictCompTime(step, curr_comp_time);
			}

		}	// Iteration Solver Loop End
		
		// HERE: Outside the iteration loop

		// ===============================================
		// ============== FINAL DATA SAVING ============== 
		// ===============================================
		
		// Summary of maximum value throughout the simulation
		_data.open("output/Parameter.dat", std::ofstream::out | std::ofstream::app);
		_data << "Maximum number of particle         :"; _data.width(12);
		_data << std::right << particle.num << "\n";

		_data << std::fixed << std::setprecision(2);
		_data << "Total computational time           :  "; _data.width(10);
		_data << std::right << cum_comp_time << " s\n";
		
		_data << std::fixed << std::setprecision(4);
		_data << "Max. global cour. number (C_max)   :  "; _data.width(10);
		_data << std::right << CourMax << "\n";
		_data << "Max. global Diff. number (phi_max) :  "; _data.width(10);
		_data << std::right << DiffMax << "\n";
		_data << "Max. global Stab. number (Re_h_max):  "; _data.width(10);
		_data << std::right << StabMax << "\n";
		_data.close();
	}
	std::cout << "\n\n#============== SIMULATION FINISHED ==============#\n\n";

	// ========== Simulation Summary ==========
	// ****************************************
	printf("+-------------- Simulation Summary ---------------+\n");
	time_t now = time(0);
	printf("Simulation end at     : %s", std::ctime(&now));
	printf("Maximum number of particle            : %8d \n", particle.num);
	printf("Total iteration number                : %8d \n", Pars::nt);
	printf("Total computational time              : %8.2f h \n", cum_comp_time/3600.0);
	printf("Average iteration computing time      : %8.2f s \n", cum_comp_time/Pars::nt);
	printf("Max. global Courant number (C_max)    : %8.2f \n", CourMax);
	printf("Max. global Diff. number (phi_max)    : %8.2f \n", DiffMax);
	printf("Max. global Stab. number (Re_h_max)   : %8.2f \n", StabMax);
	printf("+-------------------------------------------------+\n");
	std::cout << "\n#=================================================#\n\n";

	// // ADDITIONAL CODE to Show the argc and argv of main input parameter
	// std::cout << "\nThere are " << argc << " arguments:" << std::endl;
	// for (size_t i = 0; i < argc; i ++){
	// 	std::cout << argv[i] << std::endl;
	// }

	return 0;
}

/* PROGRAMMING & SIMULATION NOTEs:
   +) Makefile directory -> c:\program files\gnuwin32\bin\make.exe
   +) To avoid memory leak, please do an initalization for all struct members (e.g, set to be 0)
   +) Plese put attention to variable type definition
*/