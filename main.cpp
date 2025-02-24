/* ##########################################
 *  VORTEX PARTICLE METHOD on LSMPS disc.
 *  Sec		: main.cpp
 *
 *  Created by Angga'18 on 16/10/23. -> It takes 2 months but still need more refinement :(
 *  Co.Auth	: Adhika'15, Ical'17
 *	
 *	FOR ANY KIND OF COPY OR REDUPLICATION
 *	DO NOT DELETE THIS LABEL !!!
 * ##########################################*/

/* ---------------- PROGRAM DESCRIPTION ----------------
    > The main program of VPM-LSMPS fluid solver program
    > All parameter used in running the program is stored
      in global.cpp
    > Detailed calculation process is given later on
    > This program contains:
      - FMM accelerated biot-savart using adaptive tree cell
      - LSMPS properties interpolation
      - Brinkmann penalization method
      - Adaptive Multiresolution Particle distribution
*/

#define TESTING_POISSON false
#define TESTING_LSMPS false

// #pragma region INCLUDE_PACKAGE
	// Parameter and Utilities
	// ***********************
	#include "global.hpp"						// Simulation setting and parameter
	#include "Utils.hpp"						// Simulation data storage
	#include "src/stability/stability.hpp"		// Simulation stability criteria evaluation
	#include "src/grid_block/gridNode.hpp"		// Particle data grouping

	// Subroutine Packages
	// *******************
	#include "src/geometry/geometry.hpp"				// Body generation
	#include "src/initialization/initialization.hpp"	// Simulation domain initialization
	#include "src/remeshing/remeshing.hpp"				// A redistribution package
	#include "src/adaptation/adaptation.hpp"			// An adaptive particle distribution package
	#include "src/save_data/save_data.hpp"			    // Data write and saving package

	// Additional for LSMPS check
	// #include "src/grid_block/"			    // Data write and saving package
	#include "src/grid_block/generateGrid.hpp"
	#include "src/grid_block/gridNodeNgh.hpp"

	// Physics Calculation Method (*self explained)
	// **************************
	#include "src/advection/advection.hpp"
	#include "src/diffusion/diffusion.hpp"
	#include "src/stretching/stretching.hpp"
	#include "src/velocity_calculation/velocity_calc.hpp"
	#include "src/pressure_poisson/pressure_poisson.hpp"
	#include "src/penalization/penalization.hpp"
	#include "src/force_calculation/force_calc.hpp"
// #pragma endregion

void simUtil::timingCount(double duration, std::string nameTime, int end){
	std::ofstream _write;
	if (end == -2){
		_write.open(nameTime);
		_write << "it,intNgh,ngh,LSMPS\n";
		_write.close();
	}
	else if (end == -1){
		#if (DIM == 2)
			_write.open(nameTime);
			_write << "it,stab,pen,adv,dif,adp,rmsh,fmm\n";
			_write.close();
		#elif (DIM == 3)
			_write.open(nameTime);
			_write << "it,stab,pen,adv,difstr,adp,rmsh,fmm\n";
			_write.close();
		#endif
		// #if (DIM == 2)
		// 	_write.open(nameTime);
		// 	_write << "it,stab,adv,dif,adp,rmsh,fmm\n";
		// 	_write.close();
		// #elif (DIM == 3)
		// 	_write.open(nameTime);
		// 	_write << "it,stab,adv,difstr,adp,rmsh,fmm\n";
		// 	_write.close();
		// #endif
	}else{
		_write.open(nameTime, std::ofstream::out | std::ofstream::app);
		_write << duration;
		if (end==0){
			_write << ",";
		}else if (end==1){
			_write << "\n";
		}
		_write.close();
	}
	return;
}

#if (DATA_INTERPOLATION == 0)
#if (!TESTING_LSMPS)
// =====================================================
// +--------------- Program Starts Here ---------------+
// =====================================================
int main(int argc, char const *argv[])
{
	// ==================== Initial Definition Region ====================
	// *******************************************************************
	// #pragma region DATA_STORAGE
		// #if (N_BODY > 0) 
			std::vector<Body> bodyList(N_BODY);  // Obstacle solid object data list 
		// #else
		// 	std::vector<Body> bodyList(1);       // Obstacle solid object data list [Direct declaration, only works as helper for no body simulation]
		// #endif
		Particle particle;                   // Simulation physical data storage
		GridNode nodeGridList;               // Node data structure (Particle Related Data)
	// #pragma endregion

	// #pragma region SUBROUTINE_INSTANCES
		// *Initialization tools
		geometry geom_tool;                  // Geometry generation
		initialization initialization_tool;  // Particle distribution
		remeshing remesh_tool;               // Particle redistribution and neighbor search (Save the Eulerian basis particle groupping)
		
		// *Solver tools
		penalization penalization_tool;      // Penalization calculation
		VelocityCalc velocity_tool;      	 // Biot Savart velocity solver [accelerated by FMM]
		advection advection_tool;            // Advection calculation
		diffusion diffusion_tool;            // Diffusion calculation
		stretching stretching_tool;          // Stretching calculation
		pressure_poisson pressure_tool;	     // Poisson solver of pressure [Currently not used]
		force_calculation force_tool;        // Calculate and save force data
		
		// *Simulation tools
		simUtil utility;                     // Simulation utilities
		stability stab_eval_tool;			 // Stability evaluator
		save_data save_manager;              // Data writing manager
		std::ofstream _write;                // Data writer
		std::string nameTime = "output/Timing.csv";
		std::string nameTimeRed = "output/TimingRed.csv";
	// #pragma endregion


	// Solver computational time manager
	auto t_tick = std::chrono::system_clock::now();	// Using chrono
	std::chrono::duration<double> time_dur = (std::chrono::system_clock::now() - t_tick);
	double duration = time_dur.count();				// Current iteration computational time

	// Write the timing manager header
	if (true) utility.timingCount(0.0, nameTime, -1);
	if (true) utility.timingCount(0.0, nameTimeRed, -2);

	// #pragma region INTERNAL_VARIABLE
		double curr_comp_time = 0.0e0;       // Total computational time at each iteration
		double cum_comp_time  = 0.0e0;       // Cumulative of computational time
		int nt_start = 0;                    // The starting iteration step
			if      (Pars::opt_start_state == 0) nt_start = 0;
			else if (Pars::opt_start_state == 1) nt_start = Pars::resume_step + 1;
	// #pragma endregion

	// #pragma region SIMULATION_STABILITY_INDICATOR
		double CourMax = Pars::Courant;     // Global maximum courant number throughout simulation
		double DiffMax = Pars::Diffusion;   // Global maximum diffusion number throughout simulation
		double VortMax = 0.0;   			// Global maximum vorticity criteria throughout simulation (mesh criteria) [Ploumhans]
		double StrchMax = 0.0;   			// Global maximum stretching criteria throughout simulation (lagrangian time stability criteria)
		
		// Collection of stability number
		std::vector<double*> stabNumList;
		stabNumList.push_back(&CourMax);
		stabNumList.push_back(&DiffMax);
		stabNumList.push_back(&VortMax);
		stabNumList.push_back(&StrchMax);
	// #pragma endregion

	// [CONSOLE LOG] Display simulation parameter to console
	save_manager.summary_log();

	// [DATA  WRITE] Save the simulation parameter to file
	save_manager.write_summary_log();
	
	// Pre-processing prompt display
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RED);
	printf("+---------------- %sPRE PROCESSING%s -----------------+\n", FONT_TORQUOISE, FONT_RED);
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RESET);
	

	// =========================== Body Region ===========================
	// *******************************************************************
	// Generate all body data
	geom_tool.generateBody(bodyList);
	
	// Save all body data
	for (size_t i = 0; i < /*bodyList.size()*/N_BODY; i++)
	save_manager.save_body_state(bodyList[i], std::to_string(i+1), 0);


	// ========================= Particle Region =========================
	// *******************************************************************
	// Initial particle distribution generation
	initialization_tool.initialize_particle(particle, bodyList, nodeGridList);

	// Initialization of particle vorticity distribution
	initialization_tool.initialize_vorticity(particle, nodeGridList, bodyList);
	remesh_tool.set_neighbor(particle, nodeGridList);

	// Trim the particle data for (Lagrangian particle)
	if (Pars::flag_compact_domain && Pars::opt_start_state == 0) utility.trimParticleDomain(particle);

	// Set the initial particle data saved file name
	std::string saveNameInit;
		if 		(Pars::opt_start_state == 0) saveNameInit = "initial";	// Start at initial time
		else if (Pars::opt_start_state == 1) saveNameInit = "resume";	// Resume simulation (from the given iteration)
	
	// Save the particle data
	save_manager.save_par_state(particle, saveNameInit, 0);

	// Save the grid node data
	save_manager.save_grid_node_state(nodeGridList, saveNameInit, 1);


	// ========================= Physical Region =========================
	// *******************************************************************
	
	/* Still no code lies here ... */
	
	// Initialization the chi data
	// penalization_tool.initialize();			// [Create a chi calculation] [Already done in initialization sub.]
	// ** Vibrational Parameter Definition 	[FURTHER WORK, see most bottom part, outside the main block]
	
	
	// ===================== Pre Processing Summary ======================
	// *******************************************************************
	// Calculate the number of body node
	int sumBodyNode = 0;
	for (int i = 0; i < N_BODY; i++) sumBodyNode += bodyList[i].n_node;
	
	// Display the summary to the console
	printf("\n+------------ Initialization Summary -------------+\n");
	if (sumBodyNode == 0)	
	printf("No body simulation <!>\n");
	else					
	printf("Count of total body node                : %8d \n", sumBodyNode);
	printf("Count of particle                       : %8d \n", particle.num);
	if (Pars::opt_init_particle == 5)
	printf("Count of grid node                      : %8d \n", (int)nodeGridList.nodeMap.size());
	printf("Total iteration number                  : %8d \n", Pars::max_iter);
	printf("+-------------------------------------------------+\n\n");
	
	// Update the summary data file
	_write.open(save_manager.get_log_directory(), std::ofstream::out | std::ofstream::app);
	_write << "Number of initialized particle     :" << w12 << wR << particle.num << "\n";
	_write.close();

	// exit(0);

	// ====================== Simulation Run Prompt ======================
	// *******************************************************************
	// Simulation command prompt
	std::cout << FONT_BLUE  << "<!> The initialization is completed!\n"
			  << FONT_RESET << "<!> Do you want to run the simulation? ("
			  << FONT_GREEN << "yes" << FONT_RESET << "/" 
			  << FONT_RED   << "no"  << FONT_RESET << ")\n" 
			  << "    Type here: ";
	
	// Prompt to continue run the simulation
	// Give input to the prompt
	bool _run = false; std::string _cmd; std::cin >> _cmd;
	// Set the simulation command
	if (_cmd == "yes" || _cmd == "Yes" || _cmd == "Y" || _cmd == "y") _run = true;
	else _run = false;

	// =====================================================
	// =============== SIMULATION ITERATION ================
	// =====================================================
	if (_run == true)
	{
		// Solving header prompt
		printf("\n%s#=================================================#%s", FONT_RED, FONT_RED);
		printf("\n+-------------------- %sSOLVING%s --------------------+", FONT_TORQUOISE, FONT_RED);
		printf("\n%s#=================================================#%s", FONT_RED, FONT_RESET);

		// Particle data counter throughout simulation
		int minParticleNum = particle.num;		// Minimum number of particle throughout simulation
		int maxParticleNum = particle.num;		// Maximum number of particle throughout simulation
		
		// A standard test for
		bool peturbation_trigger = Pars::flag_peturbation; 
		bool velocity_L2N = false;
		bool vorticity_L2N = false;

		// // Data for Adam Bashforth
		// // The temporary data for the previous data
		// Particle parPrev;					// The previous particle data (to save the velocity)
		// std::vector<double> dVortXdt;		// The previous time derivation of X vorticity
		// std::vector<double> dVortYdt;		// The previous time derivation of Y vorticity
		// std::vector<double> dVortZdt;		// The previous time derivation of Z vorticity
		Particle parPrev; // The data collection at iteration n-1

		// ========================= Iteration Loop ==========================
		// *******************************************************************
		for (int step = nt_start; step < Pars::max_iter; step++){
			// ********************* ITERATION INITIALs **********************
			// Print current iteration header
			utility.printHeader(step);
			particle.currStep = step;		// Update the current step into the particle data

			if (true){
				_write.open(nameTime, std::ofstream::out | std::ofstream::app);
				_write << step << ",";
				_write.close();
				_write.open(nameTimeRed, std::ofstream::out | std::ofstream::app);
				_write << step << ",";
				_write.close();
			}

			t_tick = std::chrono::system_clock::now();	// Using chrono
			// Printing the stability criteria: courant (C), diffusion (Phi), vortex mesh reynold (Re_h), and lagrangian stretching criteria number
			if (Pars::flag_disp_stability || Pars::flag_save_stability){
				stab_eval_tool.stabilityEval(particle, stabNumList, step);
			}
			time_dur = (std::chrono::system_clock::now() - t_tick);
			duration = time_dur.count();				// Current iteration computational time
			if (true) utility.timingCount(duration, nameTime, 0);

			
			// Solver computational time manager
			auto t_start = std::chrono::system_clock::now();	// Using chrono

			// **ADDITIONAL ENSTROPHY and GLOBAL KINETIC CALCULATION
			if (Pars::flag_estrophy_calc) utility.calculate_save_enstrophy(particle, step);
			if (Pars::flag_kinetic_calc) utility.calculate_save_kinetic(particle, step);
			// utility.calculate_save_dissipation(particle,step);

			// Calculating the error
			// std::cout << "Start L2N!!\n";
			if (velocity_L2N) utility.velocity_L2Norm(particle, step, 2, 0);
			// std::cout << "Done L2N Vel!!\n";
			if (vorticity_L2N) utility.vorticity_L2Norm(particle, step, 2);
			// std::cout << "Done L2N Vort!!\n";


			// ******************** PHYSICAL CALCULATION *********************
			{
				/*
				Sequence of solving computation:
				1. Particle Redistribution  [Structured particle distribution]
				2. Velocity calculation     [Structured particle distribution]
				   -> Calculating rotational velocity
				   -> Calculating total velocity by helmholtz decomposition
				   -> Saving particle step 		[The best stage to save particle]
				3. Penalization             [Structured particle distribution]
				   -> Saving particle force 	[Must be save after penalization]
				   -> Saving particle step 		[Alternative stage]
				4. Perform Advection        [Unstructured particle distribution]
				5.1 Perform Diffusion       [Unstructured particle distribution]
				   -> Saving particle step
				5.2 Perform Stretching      [Unstructured particle distribution]
				   -> Saving particle step
				
				Try to rearrange the solver seq. 3->4->5->1->2

				Note on force calculation:
					> Penalization force calculation type must be calculated right after the penalization take place
					> The other type of force calculation still on investigation
				*/
				
				// SOLVER SEQUENCE : 3 -> 4 -> 5 -> 1 -> 2

				// // *************** ADDITIONAL STEP [1] ***************
				if (peturbation_trigger == true){
					// Set the iteration when peturbation generated
					double ptbTime = 3.0;				// The time when peturbation occured (Defined manually)
					int ptbIter = ptbTime / Pars::dt;	// The iteration step when peturbation occured
					
					// Add the peturbation into the domain
					if (step == ptbIter){
						// std::cout << "Adding a small peturbation ...\n";
						printf("%sAdding a small peturbation ...%s\n", FONT_CYAN, FONT_RESET);
						utility.addVorPertubation(particle);
						save_manager.save_par_state(particle, "peturbation", 0);	// Saving particle data	
						
						// Remove the trigger
						peturbation_trigger = false;
					}
				}
				// ***************************************************

				t_tick = std::chrono::system_clock::now();	// Using chrono
				// [2] Calculate Force for Explicit Penalization ! (Before the penalization term)
				if (Pars::opt_pen == 3)
				force_tool.force_calc(particle, bodyList, step, 2, "Exp");

				// [1] Calculate Penalization -> Stage (*,n)
				//      -> Update velocity u(*,n)
				//      -> Update vorticity vor(*,n)
				penalization_tool.get_penalization(particle, bodyList, step);
				
				// [2] Calculate Force for Implicit Penalization ! (After the penalization term)
				if (Pars::opt_pen == 1)
				force_tool.force_calc(particle, bodyList, step, 2, "Imp");

				time_dur = (std::chrono::system_clock::now() - t_tick);
				duration = time_dur.count();				// Current iteration computational time
				if (true) utility.timingCount(duration, nameTime, 0);

				// [THE SEQUENCE for Governing Equation]
				//  First order  :[1] Adv(P+C) --> [2] Diff+Str (FE)
				//  Second order :[1] Str(P) -> Adv(P) --> [2] Str(C) -> Adv(C) --> [3] Diff (FE)

				// Data saving temporary (Data saving before Predictor Corrector Stage)
				Particle parPred;
				// Assign data for the first split step calculation prediction value
				parPred.num = particle.num;
				parPred.x = particle.x;
				parPred.y = particle.y;
				parPred.s = particle.s;
				parPred.u = particle.u;
				parPred.v = particle.v;
				#if DIM == 2
				parPred.vorticity = particle.vorticity;
				#elif DIM == 3
				parPred.z = particle.z;
				parPred.w = particle.w;
				parPred.vortx = particle.vortx;
				parPred.vorty = particle.vorty;
				parPred.vortz = particle.vortz;
				#endif
				// parPred.isActive.resize(particle.num);
				// parPred.gz.resize(particle.num);

				
				t_tick = std::chrono::system_clock::now();	// Using chrono
				
				// At advection
				// ************
				advection_tool.main_advection(particle);		// Move the particle (position is change)
				advection_tool.intplt_velocity(parPred, particle, nodeGridList);  // -> Find the velocity of "particle" not the "parPred"
				advection_tool.main_advection_corr(parPred, particle);   // Move the particle (position is change)
				
				time_dur = (std::chrono::system_clock::now() - t_tick);
				duration = time_dur.count();				// Current iteration computational time
				if (true) utility.timingCount(duration, nameTime, 0);

				// // // After advection
				// // // ***************
				t_tick = std::chrono::system_clock::now();	// Using chrono
				#if (DIM == 2)
					// // FOR 2D
					diffusion_tool.main_diffusion(particle);
					diffusion_tool.diffusion_correction(parPred, particle);
				#elif (DIM == 3)
					// // FOR 3D
					stretching_tool.calc_diff_stretch(particle);
				#endif

				time_dur = (std::chrono::system_clock::now() - t_tick);
				duration = time_dur.count();				// Current iteration computational time
				if (true) utility.timingCount(duration, nameTime, 0);

				// ================================================
				// ------------- Using Adam Bashforth -------------
				// ================================================
				// Start with PC (not FE)
				if (false){
					std::cout << "THE ADAM BASHFORTH\n";
				if (step == 0 || ((step-1) % Pars::rmsh_inv) == 0 && (step != 1)){
					// Update the current particle data into the previous data
					parPrev.num = particle.num;
					parPrev.x = particle.x;
					parPrev.y = particle.y;
					parPrev.s = particle.s;
					parPrev.u = particle.u;
					parPrev.v = particle.v;
					#if DIM == 3
					parPrev.z = particle.z;
					parPrev.w = particle.w;
					parPrev.vortx = particle.vortx;
					parPrev.vorty = particle.vorty;
					parPrev.vortz = particle.vortz;
					#endif
					
					// Calculate the predictor corrector
					// Calculate the stretching
					#if (DIM == 2)
						diffusion_tool.main_diffusion(particle);
					#elif (DIM == 3)
						stretching_tool.main_stretching(particle);	// This is very simple
					#endif
					advection_tool.main_advection(particle);		// Move the particle (position is change)
					advection_tool.intplt_velocity(parPrev, particle, nodeGridList);  // -> Find the velocity of "particle" not the "parPred"
					
					#if (DIM == 2)
						diffusion_tool.diffusion_correction(parPred, particle);
					#elif (DIM == 3)
						stretching_tool.main_stretching_corr(parPrev, particle);	// This is very simple
					#endif
					advection_tool.main_advection_corr(parPrev, particle);   // Move the particle (position is change)
				}else{
					// Calculate the stretching
					#if (DIM == 2)
						diffusion_tool.diffusion_AB2(particle);
					#elif (DIM == 3)
						stretching_tool.main_stretching_AB(particle);
					#endif

					// Put the current previous data into temporary for calculation
					Particle tempPar = parPrev;

					// Update the previous particle data (For next iteration)
					parPrev.u = particle.u;
					parPrev.v = particle.v;
					#if DIM == 3
					parPrev.w = particle.w;
					#endif

					// Calculate the adam-bashforth
					advection_tool.main_advection_AB(tempPar, particle);
				}
				}

				// [4] Remesh data of vorticity
				//      -> We got new structured data
				if ((step % Pars::rmsh_inv) == 0 /*&& step != 0*/){
					remesh_tool.get_remeshing(particle, nodeGridList, step, bodyList);

					// Trim the particle data for (Lagrangian particle)
					if (Pars::flag_compact_domain) utility.trimParticleDomain(particle);
				}else if (DIM == 3){
					// Update the vorticity value (If not remeshed and in 3D simulation)
					for (int i = 0; i < particle.num; i++){
						double& vorX = particle.vortx[i];
						double& vorY = particle.vorty[i];
						double& vorZ = particle.vortz[i];
						particle.vorticity[i] = std::sqrt(vorX*vorX + vorY*vorY + vorZ*vorZ);
					}
				}

				t_tick = std::chrono::system_clock::now();	// Using chrono
				// [5] Calculate Velocity at x(n+1) and vor(n+1)
				//      -> We got new velocity [v(n+1)]
				velocity_tool.get_velocity(particle, nodeGridList, step);

				time_dur = (std::chrono::system_clock::now() - t_tick);
				duration = time_dur.count();				// Current iteration computational time
				if (true) utility.timingCount(duration, nameTime, 1);

				// [!] Calculate Force !
				// force_tool.force_calc(particle, bodyList, step, 2, "aftVel");

				// ================================================
				// --------------- Old Calculation ----------------
				// ================================================

				// // =============== SOLVER STEP [3] ===============
				// // [3] Perform the velocity penalization to induce vorticity around the body surface
				// #if (N_BODY != 0)
				// 	// [3] Perform penalization using Brinkmann: Penalize the velocity in body domain
				// 	penalization_tool.get_penalization(particle, bodyList, step);

				// 	// ================ Barrier Mark - Data Saving =================
				// 	// [!] TODO 1: Saving Data Force 
				// 	force_tool.force_calc(particle, bodyList, step, 2, "aftPen");		// Penalization mode
				// 	// force_tool.force_calc(particle, bodyList, step, 3, "aftPenMom");	// Vorticity moment mode
					
				// 	// // [!] TODO 2: Saving Particle Data
				// 	// if ((step % Pars::save_inv == 0)){
				// 	// 	// Saving particle data
				// 	// 	std::string DataName = utility.saveName(step);		// Write data file name
				// 	// 	// DataName = "aftPen_" + DataName;
				// 	// 	save_manager.save_par_state(particle, DataName, 0);	// Saving particle data
				// 	// }
				// #endif
				
				// // =============== SOLVER STEP [4] ===============
				// // [4] Convection or Advection Sub-step: Perform the particle advection
				// advection_tool.main_advection(particle);      			// [!] later: do 2nd order scheme
				
				// // =============== SOLVER STEP [5] ===============
				// // The diffusion and vortex stretching
				// #if (INVISCID_SIMULATION == 0)
				// 	// [5.1] Diffusion Sub-step: Calculate the vorticity diffusion & time integration
				// 	if (DIM == 2) diffusion_tool.main_diffusion(particle); 	// [!] later: do 2nd order scheme

				// 	// [5.2] Stretching Sub-step: Calculate the vorticity stretching & time integration
				// 	if (DIM == 3) stretching_tool.calc_diff_stretch(particle);

				// #elif (INVISCID_SIMULATION == 1)
				// 	// Stretching Sub-step: Calculate the vorticity stretching & time integration
				// 	if (DIM == 3) stretching_tool.main_stretching(particle);
				// #endif

				// // // [!] TODO 2: Saving Particle Data
				// // if ((step % Pars::save_inv == 0)){
				// // 	// Saving particle data
				// // 	std::string DataName = utility.saveName(step);		// Write data file name
				// // 	DataName = "aftStr_" + DataName;
				// // 	save_manager.save_par_state(particle, DataName, 0);	// Saving particle data
				// // }

				// // =============== SOLVER STEP [1] ===============
				// // [1] Particle redistribution: rearrange the particle distribution by interpolating vorticity
				// if ((step % Pars::rmsh_inv) == 0 || (step % Pars::adapt_inv) == 0){
				// 	// Particle redistribution (*every given iteration step)
				// 	remesh_tool.get_remeshing(particle, nodeGridList, step, bodyList);
				// }
				// // // [Additional Patch]: update the domain boundary check
				// // initialization_tool.update_domain_boundary(particle);

				// // // [!] TODO 2: Saving Particle Data
				// // if ((step % Pars::save_inv == 0)){
				// // 	// Saving particle data
				// // 	std::string DataName = utility.saveName(step);		// Write data file name
				// // 	DataName = "aftRmsh_" + DataName;
				// // 	save_manager.save_par_state(particle, DataName, 0);	// Saving particle data
				// // }

				// // =============== SOLVER STEP [2] ===============
				// // [2] Velocity calculation by biot savart: solving Rotational Velocity & Stretching
				// velocity_tool.get_velocity(particle, nodeGridList, step);

				// // [!] TODO 2: Saving Particle Data
				// if ((step % Pars::save_inv == 0)){
				// 	// Saving particle data
				// 	std::string DataName = utility.saveName(step);		// Write data file name
				// 	DataName = "aftVel_" + DataName;
				// 	save_manager.save_par_state(particle, DataName, 0);	// Saving particle data
				// }

				// // [!] DEBUG: Saving Particle Data
				// std::string DataName = utility.saveName(step);		// Write data file name
				// // DataName = "test_mres_" + DataName;
				// DataName = "hor_rec";
				// save_manager.save_par_state(particle, DataName, -1);	// Saving particle data
				
				// // Break here
				// exit(1);

			}	// End of the solver calculation

			// *********************** POST PROCESSING ***********************
			/** NOTE: The post processing contents are put right after the penalization calculation */
			// // [1] Saving Data Force 
			// if (Pars::flag_pressure_calc == true){
			// 	// Calculate the pressure
			// 	pressure_tool.get_pressure(particle);
			// }
			// force_tool.force_calc(particle, bodyList, step, 2, "force");

			// [2] Save the particle data at given step interval
			if ((step % Pars::save_inv == 0)){
				// Saving particle data
				std::string DataName = utility.saveName(step);		// Write data file name
				save_manager.save_par_state(particle, DataName, 0);	// Saving particle data	
				// save_manager.save_par_state(particle, DataName, -2);	// Saving particle data	(Analytical)
			}

			// [3] Saving Residual
			if (Pars::flag_save_residual == true)
			utility.saveResidual(particle, step);

			// [4] Update the particle count
			minParticleNum = std::min<int>(minParticleNum, particle.num);
			maxParticleNum = std::max<int>(maxParticleNum, particle.num);


			// ********************* COMPUTATIONAL TIME **********************
			// Calculate the accumulative computational time
			std::chrono::duration<double> calculation_time = (std::chrono::system_clock::now() - t_start);
			curr_comp_time = calculation_time.count();				// Current iteration computational time
			cum_comp_time += curr_comp_time;						// Accumulative computational time

			// Saving the simulation time for each iteration
			if (Pars::flag_save_sim_time){
				// Print Header
				if (step == 0 || ((Pars::opt_start_state == 1) && (step == Pars::resume_step+1))){
					_write.open("output/Simulation_Time.csv");
					_write << "iteration,sim_time,par_num,comp_time,cum_comp_time\n";
					_write.close();
				}
				
				// Print Data
				_write.open("output/Simulation_Time.csv", std::ofstream::out | std::ofstream::app);
				_write <<  "" << step
						<< "," << step * Pars::dt
						<< "," << particle.num
						<< "," << curr_comp_time
						<< "," << cum_comp_time
						<< "\n";
				_write.close();
				
			}

			// ********************** ITERATION SUMMARY **********************
			// Displaying the simulation time for each iteration
			printf("\n%s**************** Iteration Summary ****************%s", FONT_GREEN, FONT_RESET);
			printf("\n<!> Current iteration sim. time   : %12.2f s", step * Pars::dt);
			printf("\n<!> Current iteration comp. time  : %12.2f s", curr_comp_time);
			printf("\n<!> Cumulative computational time : %12.2f s\n", cum_comp_time);
			
			printf("\n<!> Particle count                : %12d", particle.num);
			printf("\n<!> Iteration to go               : %12d\n", Pars::max_iter - step);
			
			// Prediction time to finish the simulation
			if (Pars::flag_disp_pred_time == true){
				utility.predictCompTime(step, curr_comp_time);
			}
			
			// **End of the iteration
			// exit(0);
		}
		
		// HERE: Outside the iteration loop
		// printf("%s\n<!> The iteration has just finished successfully! %s\n", FONT_GREEN, FONT_RESET);


		// ======================== Final Data Saving ========================
		// *******************************************************************
		// Summary of maximum value throughout the simulation
		_write.open(save_manager.get_log_directory(), std::ofstream::out | std::ofstream::app);
		_write << "Number of particle (@t=T_end)      :" << w12 << wR << particle.num << "\n";
		_write << std::fixed << std::setprecision(2)
		       << "Total computational time           :" << w12 << wR << cum_comp_time << " s\n";
		_write << std::fixed << std::setprecision(4)
		       << "Max. global Cour. number (C_max)   :" << w12 << wR << CourMax << "\n"
		       << "Max. global Diff. number (phi_max) :" << w12 << wR << DiffMax << "\n"
		       << "Max. global Vort. number (Re_h_max):" << w12 << wR << VortMax << "\n";
		_write.close();

		
		// ======================== Simulation Summary =======================
		// *******************************************************************
		time_t now = time(0);			// Present timing

		printf("\n%s+-------------- Simulation Summary ---------------+%s\n", FONT_BLUE, FONT_RESET);
		printf("Simulation end at     : %s", std::ctime(&now));
		printf("Minimum number of particle           : %9d \n", minParticleNum);
		printf("Maximum number of particle           : %9d \n", maxParticleNum);
		printf("Total iteration number               : %9d \n", Pars::max_iter);
		
		if (cum_comp_time < 60.0){
			// If simulation runs below than 1 minute
			printf("Total computational time             : %9.2f s \n", cum_comp_time);
		}else if (cum_comp_time/60.0 < 60.0){
			// If simulation runs below than 1 hour
			int time_m = cum_comp_time/60;
			double time_s = cum_comp_time - (time_m*60);
			printf("Total computational time             : %2dm %5.2f s \n", time_m, time_s);
		}else{
			// If simulation runs longer than 1 hour
			printf("Total computational time             : %9.2f h \n", cum_comp_time/3600.0);
		}

		if (Pars::opt_start_state == 0){
			printf("Average iteration computing time     : %9.2f s \n", cum_comp_time/Pars::max_iter);
		}
		if (Pars::opt_start_state == 1){
			printf("Average iteration computing time     : %9.2f s \n", cum_comp_time/(Pars::max_iter - Pars::resume_step));
		}

		printf("Max. global Courant number (C_max)   : %9.2f \n", CourMax);
		printf("Max. global Diff. number (phi_max)   : %9.2f \n", DiffMax);
		printf("Max. global Vort. number (Re_h_max)  : %9.2f \n", VortMax);
		printf("%s+-------------------------------------------------+%s\n", FONT_BLUE, FONT_RESET);

		// Simulation is done successfully!
		printf("%s<!> The simulation is finished successfully! %s\n", FONT_GREEN, FONT_RESET);
	}
	
	else{
		// Simulation is not executed!
		printf("%s<!> The simulation is not executed! %s\n", FONT_BLUE, FONT_RESET);
	}

	return 0;
}

#else

// Additional package
#define TYPE_OF_LSMPS_DATA 3
void Util_LSMPS_gaussian_1(Particle& par);
void Util_LSMPS_gaussian_2(Particle& par);
void Util_LSMPS_gaussian_3(Particle& par);
void DATA_GENERATION(Particle& par){
	switch (TYPE_OF_LSMPS_DATA)
	{
	case 1:
		Util_LSMPS_gaussian_1(par);
		break;
	
	case 2:
		Util_LSMPS_gaussian_2(par);
		break;
	
	case 3:
		Util_LSMPS_gaussian_3(par);
		break;
	
	default:
		break;
	}
}
void Util_save_LSMPS(Particle & par, std::string name);
void Util_save_LSMPS_data(Particle & par, std::string name);

int main(int argc, char const *argv[])
{
	// ==================== Initial Definition Region ====================
	// *******************************************************************
	// #pragma region DATA_STORAGE
		#if (N_BODY > 0) 
			std::vector<Body> bodyList(N_BODY);  // Obstacle solid object data list 
		#else
			std::vector<Body> bodyList(1);       // Obstacle solid object data list [Direct declaration, only works as helper for no body simulation]
		#endif
		// Particle particle;                   // Simulation physical data storage
		// GridNode nodeGridList;               // Node data structure (Particle Related Data)
	// #pragma endregion

	// #pragma region SUBROUTINE_INSTANCES
		// *Initialization tools
		initialization initialization_tool;  // Particle distribution
		remeshing remesh_tool;               // Particle redistribution and neighbor search (Save the Eulerian basis particle groupping)
		
		// *Simulation tools
		simUtil utility;                     // Simulation utilities
		stability stab_eval_tool;			 // Stability evaluator
		save_data save_manager;              // Data writing manager
		std::ofstream _write;                // Data writer
	// #pragma endregion

	// Pre-processing prompt display
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RED);
	printf("+---------------- %sPRE PROCESSING%s -----------------+\n", FONT_TORQUOISE, FONT_RED);
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RESET);

	// Save the data for each function and tolerance type
	//   header: h, SingRe, T1, T2, T3
	// [1] Save the LSMPS L2 norm
	// [2] Save the particle number 
	// [3] Save the time for ngh
	// [4] Save the time for lsmps

	std::vector<double> tolList = {0.1,0.05,0.01};
	std::vector<int> levelMaxList = {2,3,4,5,6};
	std::vector<double> parSizeList;
	for (int i = 0; i < 5; i++) parSizeList.push_back(std::pow(10,-1.0 - (i*0.5)));

	std::string L2DataDxName = "L2NormDx.csv";
	std::string L2DataDx2Name = "L2NormDx2.csv";
	std::string L2DataDlapName = "L2NormDlap.csv";
	std::string ParNumName = "Particle.csv";
	std::string LSMPSTimeName = "LSMPSTime.csv";
	std::string NghTimeName = "NghTime.csv";
	std::vector<std::string> nameDatList = 
		{L2DataDxName, L2DataDx2Name, L2DataDlapName, 
		ParNumName, LSMPSTimeName, NghTimeName};
	std::vector<std::string> nameL2List = 
		{L2DataDxName, L2DataDx2Name, L2DataDlapName};
	// Start the saving data
	for (auto name:nameDatList){
		_write.open(name);
		_write << "h,SR,MR_T1,MR_T2,MR_T3";
		_write.close();
	}

	// Initialize the timer
	auto t_start = std::chrono::system_clock::now();	// Using chrono
	std::chrono::duration<double> time_interval;		// Timer

	// MAIN ITERATION
	for(int i = 0; i < 5; i++){
		Particle particle;                   // Simulation physical data storage
		GridNode nodeGridList;               // Node data structure (Particle Related Data)
		
		// Aliasing the particle size
		const double& parSize = parSizeList[i];
		Pars::sigma = parSize;
		Pars::max_level = levelMaxList[i];

		// Update the particle size data
		for (auto name:nameDatList){
			_write.open(name, std::ofstream::out | std::ofstream::app);
			_write << "\n" << parSize;
			_write.close();
		}

		// =======================================================================================
		// ================================== SINGLE RESOLUTION ==================================
		// =======================================================================================
		std::cout << "SINGLE RESOLUTION CALCULATION\n";
		// [1] Particle generation
		// ***********************
		initialization_tool.initialize_particle(particle, bodyList, nodeGridList);

		// Write timer
		_write.open(ParNumName, std::ofstream::out | std::ofstream::app);
		_write << "," << particle.num;
		_write.close();
		
		// [2] Neighbor evaluation
		// ***********************
		t_start = std::chrono::system_clock::now();
		remesh_tool.set_neighbor(particle, nodeGridList);
		time_interval = (std::chrono::system_clock::now() - t_start);
		double NGH_TIME = time_interval.count();

		// Write timer
		_write.open(NghTimeName, std::ofstream::out | std::ofstream::app);
		_write << "," << NGH_TIME;
		_write.close();

		// [3] Data calculation
		// ***********************
		DATA_GENERATION(particle);

		// [4] Derivative calculation
		// ***********************
		std::cout << "<-> Start calculating LSMPS Single Res\n";
		LSMPSa toolsA;
		Particle& _p = particle;

		t_start = std::chrono::system_clock::now();			// START TIMER
		toolsA.set_LSMPS(_p.x, _p.y, _p.s, _p.F_a, _p.neighbor);
		time_interval = (std::chrono::system_clock::now() - t_start);
		double LSMPS_TIME = time_interval.count();

		// Write timer
		_write.open(LSMPSTimeName, std::ofstream::out | std::ofstream::app);
		_write << "," << LSMPS_TIME;
		_write.close();

		// Retrieve the derivation data
		_p.Fx = toolsA.get_ddx();
		_p.Fy = toolsA.get_ddy();
		_p.Fx2 = toolsA.get_d2d2x();
		_p.Fy2 = toolsA.get_d2d2y();
		std::vector<double> lap(_p.num), lap_a(_p.num);
		for (int i = 0; i < _p.num; i++){
			lap[i] = _p.Fx2[i] + _p.Fy2[i];
			lap_a[i] = _p.Fx2_a[i] + _p.Fy2_a[i];
		}
		
		// L2 Norm calculation
		std::vector<double> L2Norm;
		double L2NormDx = utility.calculate_L2_Norm(_p,_p.Fx,_p.Fx_a);
		double L2NormDx2 = utility.calculate_L2_Norm(_p,_p.Fx2,_p.Fx2_a);
		double L2NormDlap = utility.calculate_L2_Norm(_p,lap,lap_a);
		L2Norm = {L2NormDx, L2NormDx2, L2NormDlap};

		// Update the particle size data
		for (int i = 0; i < 3; i ++){
			_write.open(nameL2List[i], std::ofstream::out | std::ofstream::app);
			_write << "," << L2Norm[i];
			_write.close();
		}

		// Util_save_LSMPS_data(particle,"H"+std::to_string(i+1)+"SR");
		particle.neighbor.clear();

		// =======================================================================================
		// =================================== MULTIRESOLUTION ===================================
		// =======================================================================================
		std::cout << "MULTIRESOLUTION CALCULATION\n";

		// Tolerance iteration calculation
		for(int j = 0; j < 3; j++){
			std::cout << "TOLERANCE " << j+1 << "\n";
			Particle particle_MR;
			GridNode nodeGridList_MR;
			Pars::adapt_tol = tolList[j];
			
			// // // [1] Initialization
			// // // Generate the base grid particle
			// initialization_tool.initialize_particle(particle_MR, bodyList, nodeGridList_MR);
			// DATA_GENERATION(particle_MR);

			// save_manager.save_par_state(particle_MR,"Tol"+std::to_string(j+1)+".csv",0);

			// Generate the grid from the single uniform resolution
			generateGrid genGridTools;
			genGridTools.createNode(nodeGridList_MR, particle);

			// [2] Adapt the grid, create a new multiresolution node
			adaptation adaptTools;
			neighbor nghTools;
			adaptTools.get_adaptation_LSMPS(particle, nodeGridList_MR);
			// nodeGridList_MR.saveLeafGrid(nodeGridList_MR, "Adapted");
			// Free the particle counter data
			for (auto&[_nodeID, _node] : nodeGridList_MR.nodeMap){
				_node->parList.clear();
			}

			// [3] Generate the new multiresolution particle
			// Initialize the particle that have been adapted
			Particle &par = particle_MR;
			// // Reserve the memory for particle data
			// par.x.clear(); par.u.clear(); par.vortx.clear();
			// par.y.clear(); par.v.clear(); par.vorty.clear();
			// par.z.clear(); par.w.clear(); par.vortz.clear();
			// par.gz.clear(); par.vorticity.clear();	
			// par.nodeID.clear();
			// par.level.clear();
			// par.s.clear();
			// par.chi.clear();
			// par.isActive.clear();
			// par.neighbor.clear();

			// Initialize the divisor
			int parID = 0;          // Temporary particle ID local to node
			int _parNum = Pars::intPow(nodeGridList_MR.baseParNum, DIM);  // Number of particle inside node
			int _div[DIM];          // Particle ID divisor
			int _parIndex[DIM];     // Temporary index position of particle local to node
			
			// Initialize the divisor
			_div[0] = 1;
			for (int i = 1; i < DIM; i++){
				_div[i] = _div[i-1] * nodeGridList_MR.baseParNum;
			}

			// Generate particle by loop through all nodes
			for (auto &[_ID, _node] : nodeGridList_MR.nodeMap){
				// Only do particle generation to leaf node
				if (_node->isLeaf){					
					// Calculate shared properties of particle inside the current node
					double _parSize = _node->length / nodeGridList_MR.baseParNum;

					// Generate all particle inside the curren Node
					for (int _locID = 0; _locID < _parNum; _locID++){
						// Calculate local index coordinate inside the node from local ID
						basis_loop(d) _parIndex[d] = (_locID/_div[d]) % nodeGridList_MR.baseParNum;

						// Assign the particle position
						double _x = _node->pivCoor[0] + (0.5 + _parIndex[0])*_parSize;
						par.x.push_back(_x);
						if (DIM>1) {
							double _y = _node->pivCoor[1] + (0.5 + _parIndex[1])*_parSize;
							par.y.push_back(_y);
						}
						if (DIM>2) {
							double _z = _node->pivCoor[2] + (0.5 + _parIndex[2])*_parSize;
							par.z.push_back(_z);
						}

						// Assign other data
						par.nodeID.push_back(_node->nodeID);
						par.level.push_back(_node->level);
						par.s.push_back(_parSize);

						// Insert the current particle ID into the corresponding node
						_node->parList.push_back(parID++);
					}
				}
			}

			// Update the particle number
			par.num = parID;
			{
				int &_num = par.num;
				par.gz.resize(_num); 
				par.vorticity.resize(_num);	
				par.u.resize(_num);
				par.v.resize(_num);
				par.chi.resize(_num);
				par.isActive.resize(_num);
			}

			// Re evaluate the particle
			std::cout << "\nUpdate the data set\n";
			DATA_GENERATION(particle_MR);

			// Write timer
			_write.open(ParNumName, std::ofstream::out | std::ofstream::app);
			_write << "," << particle_MR.num;
			_write.close();

			// [4] Evaluate neighbor
			t_start = std::chrono::system_clock::now();			// START TIMER
			std::cout << "\nCalculating the Neighbor\n";
			GridNodeNgh calculateToolsNgh;
			calculateToolsNgh.find_neighbor(particle_MR.neighbor, nodeGridList_MR, particle_MR);
			NGH_TIME = time_interval.count();

			// Write timer
			_write.open(NghTimeName, std::ofstream::out | std::ofstream::app);
			_write << "," << NGH_TIME;
			_write.close();

			// // Saving particle data
			// save_manager.save_par_state(particle, "After_Adapt", 0);

			// [5] Calculation on numerical derivatives
			std::cout << "<-> Start calculating LSMPS\n";
			t_start = std::chrono::system_clock::now();			// START TIMER
			toolsA.set_LSMPS(par.x, par.y, par.s, par.F_a, par.neighbor);
			time_interval = (std::chrono::system_clock::now() - t_start);
			double LSMPS_TIME = time_interval.count();

			// Write timer
			_write.open(LSMPSTimeName, std::ofstream::out | std::ofstream::app);
			_write << "," << LSMPS_TIME;
			_write.close();

			// Retrieve the derivation data
			par.Fx = toolsA.get_ddx();
			par.Fy = toolsA.get_ddy();
			par.Fx2 = toolsA.get_d2d2x();
			par.Fy2 = toolsA.get_d2d2y();
			std::vector<double> lap(par.num), lap_a(par.num);
			for (int i = 0; i < par.num; i++){
				lap[i] = par.Fx2[i] + par.Fy2[i];
				lap_a[i] = par.Fx2_a[i] + par.Fy2_a[i];
			}
			
			// L2 Norm calculation
			std::vector<double> L2Norm;
			double L2NormDx = utility.calculate_L2_Norm(par,par.Fx,par.Fx_a);
			double L2NormDx2 = utility.calculate_L2_Norm(par,par.Fx2,par.Fx2_a);
			double L2NormDlap = utility.calculate_L2_Norm(par,lap,lap_a);
			L2Norm = {L2NormDx, L2NormDx2, L2NormDlap};

			// Update the particle size data
			for (int i = 0; i < 3; i ++){
				_write.open(nameL2List[i], std::ofstream::out | std::ofstream::app);
				_write << "," << L2Norm[i];
				_write.close();
			}

			Util_save_LSMPS_data(particle_MR,"H"+std::to_string(i+1)+"T"+std::to_string(j+1));
		}
		// exit(0);
	}

// 	// How to add something into the data.csv
// 	_write.open("NAME", std::ofstream::out | std::ofstream::app);
// 	_write << "Number of particle (@t=T_end)      :" << w12 << wR << particle.num << "\n";
// 	_write.close();
	
// 	// // LSMPS TEST
// 	// // ==================================================================
// 	// bool TYPE_FULL_LSMPS = false;
// 	// // std::string segmentName = "SR_H5";	// Single resoltuion

// 	// std::string segmentName = "H3_T2";	// Multi resolution
// 	// if (TYPE_FULL_LSMPS){
// 	// 	// Generate the particle neighbor ID data list
// 	// 	remesh_tool.set_neighbor(particle, nodeGridList);
// 	// }
	
// 	// MESSAGE_LOG << "STARTING THE LSMPS CALCULATION\n";
// 	// // Basic tools and parameter
// 	// // auto t_start = std::chrono::system_clock::now();	// Using chrono
// 	// // std::chrono::duration<double> time_interval;		// Timer
// 	// LSMPSa toolsA;
// 	// Particle& par = particle;

// 	// // Saving particle data
// 	// particle.nodeID.resize(particle.num);
// 	// particle.level.resize(particle.num);
// 	// // save_manager.save_par_state(particle, "Before_Adapt", 0);
	
// 	// // Particle data generation
// 	// DATA_GENERATION(particle);	

// 	// if (TYPE_FULL_LSMPS){
// 	// 	// First calculation on spatial derivatives
// 	// t_start = std::chrono::system_clock::now();			// START TIMER
// 	// 	std::cout << "<-> Start calculating LSMPS\n";
// 	// 	toolsA.set_LSMPS(_p.x, _p.y, _p.s, _p.F_a, _p.neighbor);
// 	// 	_p.Fx = toolsA.get_ddx();
// 	// 	_p.Fy = toolsA.get_ddy();
// 	// 	_p.Fx2 = toolsA.get_d2d2x();
// 	// 	_p.Fy2 = toolsA.get_d2d2y();
// 	// time_interval = (std::chrono::system_clock::now() - t_start);		// TIMER END
// 	// double TIME_LSMPS_FULL = time_interval.count();				// Current iteration computational time

// 	// 	// Data retrieval 1
// 	// 	Util_save_LSMPS(particle, "LN_Data_" + segmentName);

// 	// 	// Print Header
// 	// 	_write.open("output/Data_summary_"  + segmentName + ".dat");
// 	// 	// Print Data
// 	// 	_write << "Data summary\n";
// 	// 	if (TYPE_OF_LSMPS_DATA == 1){
// 	// 		_write << " > Function 1: Gaussian\n";
// 	// 	}else if (TYPE_OF_LSMPS_DATA == 2){
// 	// 		_write << " > Function 2: Gaussian 2 Side\n";
// 	// 	}else if (TYPE_OF_LSMPS_DATA == 3){
// 	// 		_write << " > Function 3: Gaussian 4 Side\n";
// 	// 	}
// 	// 	_write << " > Particle Size: " << Pars::sigma << "\n"
// 	// 		<< "\n";
// 	// 	_write << "Time summary data\n"
// 	// 		<< " > LSMPS calculation : " << TIME_LSMPS_FULL << "\n"
// 	// 		<< "\n";
// 	// 	_write.close();
// 	// 	exit(0);
// 	// }
	
// 	// Evaluate AMR based on
// 	// Adaptation of --- Generate particle
// t_start = std::chrono::system_clock::now();			// START TIMER
// 	generateGrid genGridTools;
// 	genGridTools.createNode(nodeGridList, particle);
// 	// nodeGridList.saveLeafGrid(nodeGridList, "Created");

// time_interval = (std::chrono::system_clock::now() - t_start);		// TIMER END
// double TIME_GRID = time_interval.count();				// Current iteration computational time

// 	// Perform adaptation
// t_start = std::chrono::system_clock::now();			// START TIMER
// 	adaptation adaptTools;
// 	neighbor nghTools;
// 	adaptTools.get_adaptation_LSMPS(particle, nodeGridList);
// 	// nodeGridList.saveLeafGrid(nodeGridList, "Adapted");
// 	// Free the particle counter data
// 	for (auto&[_nodeID, _node] : nodeGridList.nodeMap){
// 		_node->parList.clear();
// 	}
// time_interval = (std::chrono::system_clock::now() - t_start);		// TIMER END
// double TIME_ADAPT = time_interval.count();				// Current iteration computational time
// 	// exit(0);

// 	// Initialize the particle that have been adapted
// t_start = std::chrono::system_clock::now();			// START TIMER
// 	Particle &par = particle;
// 	// Reserve the memory for particle data
//     par.x.clear(); par.u.clear(); par.vortx.clear();
//     par.y.clear(); par.v.clear(); par.vorty.clear();
//     par.z.clear(); par.w.clear(); par.vortz.clear();
// 	par.gz.clear(); par.vorticity.clear();	
//     par.nodeID.clear();
//     par.level.clear();
//     par.s.clear();
// 	par.chi.clear();
// 	par.isActive.clear();
// 	par.neighbor.clear();

//     // Initialize the divisor
//     int parID = 0;          // Temporary particle ID local to node
//     int _parNum = Pars::intPow(nodeGridList.baseParNum, DIM);  // Number of particle inside node
//     int _div[DIM];          // Particle ID divisor
//     int _parIndex[DIM];     // Temporary index position of particle local to node
    
//     // Initialize the divisor
//     _div[0] = 1;
//     for (int i = 1; i < DIM; i++){
//         _div[i] = _div[i-1] * nodeGridList.baseParNum;
//     }

//     // Generate particle by loop through all nodes
//     for (auto &[_ID, _node] : nodeGridList.nodeMap){
//         // Only do particle generation to leaf node
//         if (_node->isLeaf){
//             #if TRIM_DOMAIN_FLAG == 1
//             #else
//                 // Calculate shared properties of particle inside the current node
//                 double _parSize = _node->length / nodeGridList.baseParNum;

//                 // Generate all particle inside the curren Node
//                 for (int _locID = 0; _locID < _parNum; _locID++){
//                     // Calculate local index coordinate inside the node from local ID
//                     basis_loop(d) _parIndex[d] = (_locID/_div[d]) % nodeGridList.baseParNum;

//                     // Assign the particle position
//                     double _x = _node->pivCoor[0] + (0.5 + _parIndex[0])*_parSize;
//                     par.x.push_back(_x);
//                     if (DIM>1) {
//                         double _y = _node->pivCoor[1] + (0.5 + _parIndex[1])*_parSize;
//                         par.y.push_back(_y);
//                     }
//                     if (DIM>2) {
//                         double _z = _node->pivCoor[2] + (0.5 + _parIndex[2])*_parSize;
//                         par.z.push_back(_z);
//                     }

//                     // Assign other data
//                     par.nodeID.push_back(_node->nodeID);
//                     par.level.push_back(_node->level);
//                     par.s.push_back(_parSize);

//                     // Insert the current particle ID into the corresponding node
//                     _node->parList.push_back(parID++);
//                 }
//             #endif
//         }
//     }

//     // Update the particle number
//     par.num = parID;
// 	{
// 		int &_num = par.num;
// 		par.gz.resize(_num); 
// 		par.vorticity.resize(_num);	
// 		par.u.resize(_num);
// 		par.v.resize(_num);
// 		par.chi.resize(_num);
// 		par.isActive.resize(_num);
// 	}
	
// time_interval = (std::chrono::system_clock::now() - t_start);		// TIMER END
// double TIME_PAR_INIT = time_interval.count();				// Current iteration computational time

// 	// Calculate the neighbor
// t_start = std::chrono::system_clock::now();			// START TIMER
// 	std::cout << "\nCalculating the Neighbor\n";
// 	GridNodeNgh calculateToolsNgh;
// 	calculateToolsNgh.find_neighbor(particle.neighbor, nodeGridList, particle);
// time_interval = (std::chrono::system_clock::now() - t_start);		// TIMER END
// double TIME_NGH = time_interval.count();				// Current iteration computational time

// 	// // Saving particle data
// 	// save_manager.save_par_state(particle, "After_Adapt", 0);
    
// 	// GridNodeNgh nghTools;
// 	// Generate the node grid
	
// 	// Update the node grid following the adaptive multiresolution
// 	// Calculate the neighbor

// 	// Re evaluate the particle
// 	std::cout << "\nUpdate the data set\n";
// 	DATA_GENERATION(particle);

// 	// Final calculation on spatial derivatives
// 	// First calculation on spatial derivatives
// 	std::cout << "<-> Start calculating LSMPS\n";
// t_start = std::chrono::system_clock::now();			// START TIMER
// 	toolsA.set_LSMPS(_p.x, _p.y, _p.s, _p.F_a, _p.neighbor);
// 	_p.Fx = toolsA.get_ddx();	
// 	_p.Fy = toolsA.get_ddy();
// 	_p.Fx2 = toolsA.get_d2d2x();
// 	_p.Fy2 = toolsA.get_d2d2y();

// time_interval = (std::chrono::system_clock::now() - t_start);		// TIMER END
// double TIME_LSMPS = time_interval.count();				// Current iteration computational time

// 	// // Solver computational time manager
// 	// auto t_start = std::chrono::system_clock::now();	// Using chrono
// 	// std::chrono::duration<double> time_interval = (std::chrono::system_clock::now() - t_start);
// 	// curr_comp_time = time_interval.count();				// Current iteration computational time
	
// 	// Data retrieval
// 	// Util_save_LSMPS_data(particle, "After_Adapt");
// 	Util_save_LSMPS(particle, "LN_Data_MR_" + segmentName);

// 	// Print Header
// 	_write.open("output/Data_summary_MR_" + segmentName + ".dat");
// 	// Print Data
// 	_write << "Data summary\n";
// 	if (TYPE_OF_LSMPS_DATA == 1){
// 		_write << " > Function 1: Gaussian\n";
// 	}else if (TYPE_OF_LSMPS_DATA == 2){
// 		_write << " > Function 2: Gaussian 2 Side\n";
// 	}else if (TYPE_OF_LSMPS_DATA == 3){
// 		_write << " > Function 3: Gaussian 4 Side\n";
// 	}
// 	_write << " > Particle size     : " << Pars::sigma     << "\n"
// 		   << " > Max level         : " << Pars::max_level << "\n"
// 		   << " > Adaptation tol.   : " << Pars::adapt_tol << "\n"
// 		   << "\n";
// 	_write << "Time summary data\n"
// 		   << " > Grid generation   : " << TIME_GRID     << "\n"
// 		   << " > Adaptation grid   : " << TIME_ADAPT    << "\n"
// 		   << " > Particle initial  : " << TIME_PAR_INIT << "\n"
// 		   << " > Neighbor search   : " << TIME_NGH      << "\n"
// 		   << " > LSMPS calculation : " << TIME_LSMPS    << "\n"
// 		   << "\n";
// 	_write.close();
	
// 	exit(0);
// 	// ==================================================================
	return 0;
}
#endif

#elif (DATA_INTERPOLATION == 1)
#include "src/grid_block/generateGrid.hpp"
// =====================================================
// +--------------- Interpolation Data ----------------+
// =====================================================
int main(int argc, char const *argv[])
{
	/** 
	 *  @brief	A main subroutine for particle data interpolation. This process belongs to 
	 *  post processing subroutines of multiresolution data in 3D space. In order to 
	 *  efficiently render the data into paraview, a structured single resolution data is 
	 *  necessary. This subroutine will process all saved particle state data from the
	 *  previously computed simulation. The interpolated data are labeled by "SRID".
	 * 
	 *  NOTE:
	 * 	 > This subroutine only belong to 3D space simulation, since 2D simulation is not necessary.
	 *   > Subroutine also provided with Q and λ2 calculated data.
	 *   > The collected data is 
	 * 		(1) Coordinate[3]
	 * 		(2) Velocity[3]
	 * 		(3) Vorticity[3]
	 * 		(4) Q[1]
	 * 		(5) λ2[1]
	 *     of total 11 variable set
	 * 
	 *  PREPARE:
	 *   > Set the dimension in 3D (Still need to be adjusted)
	 *   > Define the target domain size and interpolated particle spacing (see global.cpp)
	 *   > 
	*/

	// Main interpolation parameter (can be adjusted manually)
	// const int fin_iter = Pars::max_iter;            // The final iteration for interpolation
	// const int save_int = Pars::save_inv;			// Iteration interval for saving data
	// const int data_num = 1 + fin_iter/save_int;     // The number of particle state data (including the zeros)

	const int fin_iter = 2000; // 24000;            // The final iteration for interpolation
	const int save_int = 40;   // 20;			// Iteration interval for saving data
	const int data_num = 1 + fin_iter/save_int;     // The number of particle state data (including the zeros)

	// #pragma region SUBROUTINE_INSTANCES
		initialization initialization_tool;  // Particle distribution
		remeshing remesh_tool;               // Particle redistribution
		simUtil utility;                     // Simulation utilities
		save_data save_manager;              // Data writing manager
	// #pragma endregion

	// Pre-processing prompt display
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RED);
	printf("+-------------------- %sSUMMARY%s --------------------+\n", FONT_TORQUOISE, FONT_RED);
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RESET);

	int xCnt = std::ceil(Pars::lxdomInt/Pars::sigmaInt);
	int yCnt = std::ceil(Pars::lydomInt/Pars::sigmaInt);
	int zCnt = std::ceil(Pars::lzdomInt/Pars::sigmaInt);

	// Initial flow parameters summary data
	printf("\n+--------- Interpolation Parameters Data ---------+\n");
	printf(" Domain x length                       : %7.2f m\n", Pars::lxdomInt);
	printf(" Domain y length                       : %7.2f m\n", Pars::lydomInt);
	printf(" Domain z length                       : %7.2f m\n", Pars::lzdomInt);
	printf(" Interpolated core size                : %7.2f m\n", Pars::sigmaInt);
	printf(" Number of particle in x basis         : %9d m\n", xCnt);
	printf(" Number of particle in y basis         : %9d m\n", yCnt);
	printf(" Number of particle in z basis         : %9d m\n", zCnt);
	printf(" Total particle number                 : %9d m\n", xCnt*yCnt*zCnt);
	printf(" Final iteration step                  : %9d \n", fin_iter);
	printf(" Step interval                         : %9d \n", save_int);
	printf(" Number of sate data                   : %9d \n", data_num);
	printf("+-------------------------------------------------+\n\n");


	// ====================== Simulation Run Prompt ======================
	// *******************************************************************
	// Simulation command prompt
	std::cout << FONT_RESET << "<!> Proceed to the interpolation? ("
			  << FONT_GREEN << "yes" << FONT_RESET << "/" 
			  << FONT_RED   << "no"  << FONT_RESET << ")\n" 
			  << "    Type here: ";
	
	// Prompt to continue run the simulation
	// Give input to the prompt
	bool _run = false; std::string _cmd; std::cin >> _cmd;
	// Set the simulation command
	if (_cmd == "yes" || _cmd == "Yes" || _cmd == "Y" || _cmd == "y") _run = true;
	else _run = false;


	// ================ Iteration for data interpolation =================
	// *******************************************************************
	if (_run == true)
	{
		// Print banner
		printf("\n%s#=================================================#%s", FONT_RED, FONT_RED);
		printf("\n+----------------- %sINTERPOLATION%s -----------------+", FONT_TORQUOISE, FONT_RED);
		printf("\n%s#=================================================#%s", FONT_RED, FONT_RESET);
		
		// Computational time accumulation
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			double _time = omp_get_wtime();
		#endif
		
		
		// ======================== Interpolation Loop =======================
		// *******************************************************************
		// Start the interpolation to all particle data state
		for (int n = 0; n < data_num; n++){
			// Solver computational time manager
			auto t_start = std::chrono::system_clock::now();	// Using chrono

			// Alias to the current iteration* [Can be adjusted manually]
			const int step = n * save_int;
			
			// Print header
			utility.printHeader(step);

			// Internal variable
			Particle srcParticle;      // Source particle data container (particle from local storage)
			Particle intParticle;      // Interpolated particle data container
			GridNode nodeGridList;     // Node data structure (Related to the source particle)

			// Read the source data
			#if (DIM == 2)
				initialization_tool.read_2d_particle(srcParticle, step);
			#elif (DIM == 3)
				initialization_tool.read_3d_particle(srcParticle, step);
			#endif

			// Create grid of source data
			printf("%s\nGenerate grid block ... %s\n", FONT_CYAN, FONT_RESET);
			generateGrid grid_step;
			grid_step.createNode(nodeGridList, srcParticle);
			
			// Perform interpolation
			remesh_tool.re_arrange_distribution(intParticle, srcParticle, nodeGridList);
			
			// Save the interpolated data
			std::string stepName = utility.saveName(step);
			save_manager.save_par_interpolation(intParticle, stepName);

			// Calculate the accumulative computational time
			std::chrono::duration<double> calculation_time = (std::chrono::system_clock::now() - t_start);
			double curr_comp_time = calculation_time.count();				// Current iteration computational time

			// Prompt for display
			printf("\n%s**************** Iteration Summary ****************%s", FONT_GREEN, FONT_RESET);
			printf("\n<!> Current iteration comp. time  : %12.2f s", curr_comp_time);			
			printf("\n<!> Particle count                : %12d", intParticle.num);
			printf("\n<!> Iteration to go               : %12d\n", data_num - n);

			// Prediction time to finish the simulation
			if (Pars::flag_disp_pred_time == true){
				// Internal variable
				int est_time_d, est_time_h, est_time_m; double est_time_s;
				est_time_s = curr_comp_time * double(data_num - n - 1);
				
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
				printf("\n<!> Estimation time to finish run : %12.2f s", curr_comp_time*double(data_num - step));
				if (est_time_d == 0){
					printf("\n<!> Estimation time to finish run :    %2dh %2dm %2ds", est_time_h, est_time_m, (int)est_time_s);
				}else{
					printf("\n<!> Estimation time to finish run :  %3ddays %2dhrs", est_time_d, est_time_h);
				}

				// printf("\n");
			}
		}

		// Interpolation summary!
		// Particle generation summary time display
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
			double _time = span.count();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime() - _time;
		#endif

		
		// ======================== Simulation Summary =======================
		// *******************************************************************
		time_t now = time(0);			// Present timing

		printf("\n%s+-------------- Simulation Summary ---------------+%s\n", FONT_BLUE, FONT_RESET);
		printf("Interpolation end at  : %s", std::ctime(&now));
		if (_time < 60.0){
			// If simulation runs below than 1 minute
			printf("Total computational time             : %9.2f s \n", _time);
		}else if (_time/60.0 < 60.0){
			// If simulation runs below than 1 hour
			int time_m = _time/60;
			double time_s = _time - (time_m*60);
			printf("Total computational time             : %2dm %5.2f s \n", time_m, time_s);
		}else{
			// If simulation runs longer than 1 hour
			printf("Total computational time             : %9.2f h \n", _time/3600.0);
		}
		printf("Average iteration computing time     : %9.2f s \n", _time/data_num);
		

		// Simulation is done successfully!
		printf("%s\n<!> The post-processing is finished successfully! %s\n", FONT_GREEN, FONT_RESET);
	}
	else{
		// Interpolation is not executed!
		printf("%s\n<!> The interpolation is cancelled! %s\n", FONT_BLUE, FONT_RESET);
	}
	

	return 0;
}
#endif

/* PROGRAMMING & SIMULATION NOTEs:
   +) Makefile directory -> c:\program files\gnuwin32\bin\make.exe
   +) To avoid memory leak, please do an initalization for all struct members (e.g, set to be 0)
   +) Plese put attention to variable type definition
   +) Remember to free the memory of malloc variable
*/

/*	List of task NEED TO BE DONE
	> Generate particle		[DONE]
	> Adaptation method		[DONE] Also done the linking, but not the cell list, (further a do? The answer is not)
	> Link the remeshing 	[DONE]
	> Redistribution 		[DONE] -> CHECK LSMPS (have been adjusted with paralel)
	> Neigbor evaluation	[DONE]
	> Set a paralel computing  [ALMOST DONE] -> There are something strange things by using parallel
		>>> Syntax	: #pragma omp parallel for
		- Too high number of cores used creates an unstable sequence, break ups the simulation data
		- At some point (a parallel loop deep inside the loop), the parallel make a sequence hold so it takes longer time than usual instead
	> Check all physical subroutine
		- Velocity ? (Check the vector data size)
		- Penalization (adjust to the number of body part)
		- Advection, Diffusion, just need to check...
*/

/*
	// DEBUGING PROCESS
	// Put the variable adaptive as 
	particle.vorticity.clear();
	particle.vorticity.resize(particle.num,1e10);
	for (Body &_body : bodyList){
		geom_tool.distance_calc(particle,_body,true);
		for (int i = 0; i < particle.num; i++){
			particle.vorticity[i] = std::min<double>(particle.R[i], particle.vorticity[i]);
		}
	}
	for (int i = 0; i < particle.num; i++){
		particle.vorticity[i] = std::pow(1.0 / (1 + particle.vorticity[i]), 2.0);
	}

	// for (int i = 0; i < particle.num; i++){
	// 	particle.vorticity[i] = 0.0;
	// }

	// particle.vorticity[60] = 1000.0;

	save_manager.save_par_state(particle,"Before_Adapted", 0);

	remesh_tool.get_remeshing(particle,nodeGridList,0);
	
	// adaptation adapt_step;
	// particle.neighbor.resize(particle.num);
	// adapt_step.get_adaptation(particle,&particle,nodeGridList);
	// adapt_step.get_new_particle(particle);
	// }

	particle.u.clear();
	particle.v.clear();
	particle.chi.clear();
	particle.gz.clear();
	particle.u.resize(particle.num);
	particle.v.resize(particle.num);
	particle.chi.resize(particle.num);
	particle.gz.resize(particle.num);

	particle.vorticity.clear();
	particle.vorticity.resize(particle.num,1e10);
	for (Body &_body : bodyList){
		geom_tool.distance_calc(particle,_body, true);
		for (int i = 0; i < particle.num; i++){
			particle.vorticity[i] = std::min<double>(particle.R[i], particle.vorticity[i]);
		}
	}
	for (int i = 0; i < particle.num; i++){
		particle.vorticity[i] = std::pow(1.0 / (1 + particle.vorticity[i]), 2.0);
		particle.s[i] = Pars::sigma * Pars::intPow(2,Pars::max_level-particle.level[i]);
		// particle.vorticity[i] = 1.0;
	}
*/

/* Vibrational Parameter Definition [FURTHER WORK]
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