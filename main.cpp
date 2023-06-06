#pragma region include
	#ifndef INCLUDED_UTILS		// Already include global
	#include "Utils.hpp"
	#endif

	#ifndef INCLUDED_DATA_SAVING
	#include "src/save_data/data_saving.hpp"
	// Data stored at each or some iteration:
		// The simulation parameter (number, domain size, etc.) 	>>> make a log
		// Run time counter (at each iteration) 		    		>>> Data list per count iteration
		// Particle data (Neighbor list, chi, active particle, etc.)>>> Data list per count iteration
		// Property contour data (V, w, P, T)						>>> Data list per count iteration
		// Additional data (drag, force, moment, lift)				>>> Data list per count iteration
	#endif

	#ifndef INCLUDED_SAVE_DATA_BASE
	#include "src/save_data/base_save_data.hpp"
	#endif

	#ifndef INCLUDED_GEOMETRY
	#include "src/geometry/geometry.hpp"
	#endif

	#ifndef INCLUDED_INITIALIZATION
	#include "src/initialization/initialization.hpp"
	#endif

	#ifndef INCLUDED_REMESHING
	#include "src/remeshing/remeshing.hpp"
	#endif

	#ifndef INCLUDED_ADVECTION
	#include "src/advection/advection.hpp"
	#endif

	#ifndef INCLUDED_VELOCITY_POISSON
	#include "src/velocity_poisson/velocity_poisson.hpp"
	#endif

	#ifndef INCLUDED_DIFFUSION
	#include "src/diffusion/diffusion.hpp"
	#endif

	#ifndef INCLUDED_PENALIZATION
	#include "src/penalization/penalization.hpp"
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

	#ifndef INCLUDED_PRESSURE_POISSON
	#include "src/pressure_poisson/pressure_poisson.hpp"
	#endif


#pragma endregion

// ************* Program Starts Here ************* //
int main(int argc, char const *argv[])
{

#pragma region data_storage	        // Creating data storage object instance
	Body body;		    // The solid object
	Particle particle; 	// The simulation particle object
#pragma endregion

#pragma region instances 	        // Creating instance of subroutine
	// Initialization tools
	geometry geom_step;                  // Geometry generation
	initialization initialization_step;  // Particle distribution
	remeshing remesh_step;               // Particle redistribution and neighbor search
	
	// Solver tools
	penalization penalization_step;      // Penalization
	velocity_poisson velocity_step;      // FMM:Biot Savart velocity solver
	advection advection_step;       	 // Advection
	diffusion diffusion_step;            // Diffusion
	solverpoisson solver_poisson;        // Poisson solver
	pressure_poisson pressure_step;		 // Poisson solver of pressure
	
	// Data saving tools
	save_data save_step;                 // Save particle distribution
	base_save_data d_base_save_data;     // Save force each time interval
#pragma endregion

#pragma region summary_log	        // Displaying and saving parameter summary
	save_step.summary_log();
#pragma endregion

#pragma region internal_variables	// This is used for continuing simulation
	double curr_comp_time = 0.0e0;  	// Total computational time at each iteration
	double cum_comp_time = 0.0e0;   	// The total cumulative of simulation computational time
	int particle_max_number = 0;    	// The maximum number of particle throughout the simulation
	int nt_start = 40;	            	// The starting iteration step
	int maxDigitLen = 1 + std::floor(std::log10(Pars::nt));		// Digit of the maximum iteration
	double CourMax = Pars::Courant;	    // The global maximum courant number throughout simulation
	double DiffMax = Pars::Diffusion;	// The global maximum diffusion number throughout simulation
#pragma endregion
	
	// Pre processing prompt displaying
	printf("\n#=================================================#");
	printf("\n+---------------- PRE PROCESSING -----------------+");
	printf("\n#=================================================#\n");

	// ========== Body Region ==========
	geom_step.generateBody(body);			// Generate the body data
	save_step.save_state(body,"initial");	// Save the body data

	// ========== Particle Region ==========
	// Initial particle distribution generation
	if (Pars::opt_cont == 0){
		initialization_step.generate(body, particle);
		nt_start = -1;
	}else if(Pars::opt_cont == 1){
		initialization_step.continue_simulation(body, particle, nt_start);   // Read the last particle data
	}
	
	// Generate cell list and perform neighbor search
	remesh_step.remeshing_init(particle);
	
	particle.P.resize(particle.num,0.0);
	// Save the particle data
	if (Pars::opt_cont == 0){
		save_step.save_state(particle,"initial", 2);
	}else if (Pars::opt_cont == 1){
		save_step.save_state(particle,"continue", 2);
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
	printf("\n+------------ Initialization Summary -------------+\n");
	printf("Number of body node                     : %8d \n", body.num);
	printf("Number of particle node                 : %8d \n", particle.num);
	printf("Iteration number                        : %8d \n", Pars::nt);
	printf("+-------------------------------------------------+\n");
	
	// ========== Update Initialization Summary ==========
	std::ofstream _data;
	_data.open("output/Parameter.dat", std::ofstream::out | std::ofstream::app);
	_data << "Number of initialized particle     :"; _data.width(12);
	_data << std::right << particle.num << "\n";
	_data.close();

	// ========== Simulation Run Prompt ==========
	// Simulation command prompt
	std::cout << std::endl <<
	             "<!> The initialization is completed!\n" << 
	             "<!> Do you want to run the simulation? (yes/no)\n" <<
				 "    Type here: ";
    
	// Simulation command Input
	bool _run = false;
	std::string _cmd;
	std::cin >> _cmd;

	// Run or don't run this simulation
	if (_cmd == "yes" | _cmd == "Yes")   // RUN
	{
		_run = true;
	}
	else                                 // DON'T RUN
	{
		_run = false;
	}

	// ===============================================
	// =============== SIMULATION RUN ================
	// ===============================================
	if (_run == true)
	{
		// Solving processing header prompt
		printf("\n#=================================================#");
		printf("\n+-------------------- SOLVING --------------------+");
		printf("\n#=================================================#");		
		bool resetFMMTree = false;
		
		// *************** Iteration Loop Starts Here *************** //
		for (size_t step = nt_start + 1; step < Pars::nt; step++)
		{
			// *************** ITERATION INITIAL SETTING and DISPLAY *************** //
			// Calculate the number of iteration digit
			int iterDigitLen = 1;
			if (step != 0)iterDigitLen = 1 + std::floor(std::log10(step));

			// Printing the iteration HEADER
			{
				printf("\n\n+");
				for(int _i = 0; _i < 16 - iterDigitLen/2; _i ++){
					printf("-");
				}
				printf(" Iteration Step %d ", (int)step);
				for(int _i = 0; _i < 16 - iterDigitLen/2 - iterDigitLen%2; _i ++){
					printf("-");
				}
				printf("+\n");
			}

			// Printing the stability criteria: courant (C) and diffusion (Phi) number
			if (Pars::flag_disp_stability){
				// Initialize parameter
				double _courNum [3] = {0,0,0};  // [Current, Accumulative, Maximum]
				double _diffNum [3] = {0,0,0};  // [Current, Accumulative, Maximum]
				for (int _i = 0; _i < particle.num; _i++){
					// Calculating courant number
					_courNum[0] = std::sqrt(std::pow(particle.u[_i],2) + std::pow(particle.v[_i],2)) * Pars::dt / particle.s[_i];
					_courNum[1] += _courNum[0];
					_courNum[2] = _courNum[2] > _courNum[0] ? _courNum[2] : _courNum[0];

					// Calculating diffusion number
					_diffNum[0] = Pars::NU * Pars::dt / (particle.s[_i] * particle.s[_i]);
					_diffNum[1] += _diffNum[0];
					_diffNum[2] = _diffNum[2] > _diffNum[0] ? _diffNum[2] : _diffNum[0];
				}
				// Average value
				_courNum[0] = _courNum[1] / particle.num;
				_diffNum[0] = _diffNum[1] / particle.num;
				
				// Displaying the value
				printf("Average courant number (C_av)           : %8.4f \n", _courNum[0]);
				printf("Average diffusion number (Phi_av)       : %8.4f \n", _diffNum[0]);
				printf("Max courant number (C_max)              : %8.4f \n", _courNum[2]);
				printf("Max diffusion number (Phi_max)          : %8.4f \n", _diffNum[2]);

				// Update the global stability criteria maximum value
				if (_courNum[2] > CourMax){CourMax = _courNum[2];}
				if (_diffNum[2] > DiffMax){DiffMax = _diffNum[2];}
			}

			// Solver computational time manager
			auto t_start = std::chrono::system_clock::now();

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
				if ((Pars::opt_cont == 1) && (step == nt_start + 1)){__step = 0;}else{__step = step;}
				
				penalization_step.get_penalization(particle, body, __step);
				// Iterative penalization calculation
				if(Pars::iterative == 2){
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
				
				/*
				// DATA STATE SAVING after PENALIZATION
				if (step % Pars::step_inv == 0){
					int addDigit = maxDigitLen - iterDigitLen;
					
					// Writing data file name
					std::string DataName = "";
					DataName.append("PEN");
					for (int _spc = 0; _spc < addDigit; _spc++)
						DataName.append("0");
					DataName.append(to_string(step));
					
					// Saving particle data
					save_step.save_state(particle, DataName, 2);
				}
				*/
				
				// // =============== SOLVER STEP [2] ===============
				// // [2] Velocity calculation by poisson/biot savart: solving Rotational Velocity & Stretching
				// // int __step;
				// if (((Pars::opt_cont == 1) && (step == nt_start + 1)) || resetFMMTree){__step = 0;}else{__step = step;}
				// // std::cout << "The value of _step : " << __step << "\n";
				// velocity_step.get_velocity(particle, __step);
				// // Helmholtz decomposition
				// // u = u_rotational + u_irrotational
				// printf("<+> Adding the irrotational velocity term \n");
				// printf("    [helmholtz backward decomposition]\n");
				// for (size_t i = 0; i < particle.num; i++){
				// 	particle.u[i] += Pars::u_inf;
				// 	particle.v[i] += Pars::v_inf;
				// }

				// ================ Barrier Mark =================
				// [!] TODO: Saving Data Force 
				// Force calculation and saving [Inserted after the Penalization calculation]
				std::cout << "\nSaving force data ...\n";
				if (Pars::force_type == 1){
					std::cout << "<+> Penalization Force\n";
					double xpus[2] = {0, 0};
					d_base_save_data.force_pen(step, particle, xpus);
				}
				else if (Pars::force_type == 2){
					save_step.Force2(step,1,2,3,4,1,2,particle);
				}
				
				// =============== SOLVER STEP [4] ===============
				// [4] Convection/Advection Sub-step: Perform the particle advection
				advection_step.advection_euler(particle);      // ! later: do 2nd order scheme

				/*
				// DATA SAVING after DIFFUSION
				if (step % Pars::step_inv == 0){
					int addDigit = maxDigitLen - iterDigitLen;
					
					// Writing data file name
					std::string DataName = "";
					DataName.append("ADV");
					for (int _spc = 0; _spc < addDigit; _spc++)
						DataName.append("0");
					DataName.append(to_string(step));
					
					// Saving particle data
					save_step.save_state(particle, DataName, 2);
				}
				*/
				
				// =============== SOLVER STEP [5] ===============
				// [5] Diffusion Sub-step: Calculate the vorticity diffusion & time integration
				diffusion_step.main_diffusion(particle); // ! later: do 2nd order scheme
				
				/*
				// DATA SAVING after DIFFUSION
				if (step % Pars::step_inv == 0){
					int addDigit = maxDigitLen - iterDigitLen;
					
					// Writing data file name
					std::string DataName = "";
					DataName.append("DIF");
					for (int _spc = 0; _spc < addDigit; _spc++)
						DataName.append("0");
					DataName.append(to_string(step));
					
					// Saving particle data
					save_step.save_state(particle, DataName, 2);
				}*/
				

				// =============== SOLVER STEP [1] ===============
				// [1] Particle redistribution: rearrange the particle distribution by interpolating vorticity
				if ((step % Pars::nrmsh) == 0) // redistribute particle every given iteration step
				{
					resetFMMTree = remesh_step.get_remeshing(particle, body, step);
				}

				/*
				if (step % Pars::step_inv == 0){
					// Writing data file name
					std::string DataName = "RMSH_";
					int addDigit = maxDigitLen - iterDigitLen;
					for (int _spc = 0; _spc < addDigit; _spc++)
						DataName.append("0");
					DataName.append(to_string(step));
					
					// Saving particle data
					save_step.save_state(particle, DataName, 2);
				}
				*/

				// =============== SOLVER STEP [2] ===============
				// [2] Velocity calculation by poisson/biot savart: solving Rotational Velocity & Stretching
				// int __step;
				if (((Pars::opt_cont == 1) && (step == nt_start + 1)) || resetFMMTree){__step = 0;}else{__step = step;}
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

				// Calculating pressure
				// if (Pars::flag_pressure_calc){
				// 	pressure_step.get_pressure(particle);
				// }
				
				// if (step > 20){
				// 	pressure_step.get_pressure(particle);
				// }

			}	// End of the solver calculation

			// *************** SAVING DATA *************** //
			// Save the particle data at given step interval
			if (step % Pars::step_inv == 0){
				if (Pars::flag_pressure_calc){
					pressure_step.get_pressure(particle);
				}
				
				// Writing data file name
				std::string DataName;
				int addDigit = maxDigitLen - iterDigitLen;
				for (int _spc = 0; _spc < addDigit; _spc++)
					DataName.append("0");
				DataName.append(to_string(step));
				
				// Saving particle data
				save_step.save_state(particle, DataName, 2);
			}

			// Calculate the accumulative computational time
			std::chrono::duration<double> elapsed_time_ms = (std::chrono::system_clock::now() - t_start);
			curr_comp_time = elapsed_time_ms.count();
			cum_comp_time += curr_comp_time;

			// Saving the simulation time for each iteration
			if (Pars::flag_save_sim_time){
				if (step == 0){
					_data.open("output/Simulation_Time.csv");
					_data << "iteration,sim_time,comp_time,curr_cum_comp_time\n";
					_data << step << "," << step * Pars::dt << "," << curr_comp_time << "," << cum_comp_time << "\n";
					_data.close();
				}else{
					_data.open("output/Simulation_Time.csv", std::ofstream::out | std::ofstream::app);
					_data << step << "," << step * Pars::dt << "," << curr_comp_time << "," << cum_comp_time << "\n";
					_data.close();
				}
			}

			// Displaying the simulation time for each iteration
			printf("\n<!> Current simulation time:           %9.3f s \n", step * Pars::dt);
			printf("\n<!> Current cumulative comp. time:     %9.3f s", curr_comp_time);
			printf("\n<!> Total cumulative comp. time:       %9.3f s", cum_comp_time);

			// Prediction time to finish the simulation
			bool _predTime = true;
			if (_predTime == true){
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
			}

			// Update the maximum particle number
			particle_max_number = particle.num > particle_max_number ? particle.num : particle_max_number;

		}	// Iteration Solver Loop End
		
		// HERE: Outside the iteration loop

		// ===============================================
		// ============== FINAL DATA SAVING ============== 
		// ===============================================
		
		// Summary of maximum value throughout the simulation
		_data.open("output/Parameter.dat", std::ofstream::out | std::ofstream::app);
		_data << "Maximum number of particle         :"; _data.width(12);
		_data << std::right << particle_max_number << "\n";
		_data << std::fixed << std::setprecision(2);
		_data << "Total computational time           :  "; _data.width(10);
		_data << std::right << cum_comp_time << " s\n";
		_data << std::fixed << std::setprecision(4);
		_data << "Max. global cour. number (C_max)   :  "; _data.width(10);
		_data << std::right << CourMax << "\n";
		_data << "Max. global Diff. number (phi_max) :  "; _data.width(10);
		_data << std::right << DiffMax << "\n";
		_data.close();
	}

	std::cout << "\n\n#============== SIMULATION FINISHED ==============#\n\n";
	
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