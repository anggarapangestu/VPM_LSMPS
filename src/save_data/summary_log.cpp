#include "data_saving.hpp"

void save_data::summary_log(){
	// Simulation Log Initial Parameter Summary Data
	std::cout << "#=================================================#\n";
	std::cout << "+---------------- SIMULATION LOG -----------------+\n";
	std::cout << "#=================================================#\n";
	std::cout << "<!> Computation Message Log\n";
	
	// Simulation parameter summary data
	printf("\n+---------- Simulation Parameters Data -----------+\n");
	printf("Body option                           :   type %d \n", Pars::opt_body);
	printf("Initialization option                 :   type %d \n", Pars::opt_init);
	printf("Neighbor search option                :   type %d \n", Pars::opt_neighbor);
	printf("Penalization option                   :   type %d \n", Pars::opt_pen);
	printf("Maximum resolution level              :        %d \n", Pars::max_level);
	printf("Core size                             : %8.4f m\n", Pars::sigma);
	printf("Time step                             : %8.4f s\n", Pars::dt);
	printf("Total simulation time                 : %8.2f s\n", Pars::sim_time);
	printf("Total iteration step                  : %8d\n", Pars::nt);
	printf("+-------------------------------------------------+\n");
	
	// Initial flow parameters summary data
	printf("\n+------------- Flow Parameters Data --------------+\n");
	printf("Reynolds number (RE)           : %10.2f [-]\n", Pars::RE);
	printf("Freestream velocity (U)        : %10.2f m/s\n", Pars::U_inf);
	printf("Fluid density (rho)            : %10.2f kg/m^3\n", Pars::RHO);
	printf("Fluid viscosity (nu)           : %10f m^2/s\n", Pars::NU);
	printf("+-------------------------------------------------+\n");

	// Additional calculation data
	printf("\n+---------------- Additional Data ----------------+\n");
	printf("Courant number (C < 1.0)            : %12f \n", Pars::Courant);
	printf("Diffusion number (Phi < 0.5)        : %12f \n", Pars::Diffusion);
	printf("Stability Criteria (Re_h < 1.0)     : %12f \n", 100 * Pars::sigma*Pars::sigma/Pars::NU);
	printf("Turbulent scaling                   : %12f \n", std::ceil(1.0e0 / Pars::Courant));
	printf("Max time step (Re_h criteria)       : %12f \n", 0.25*Pars::sigma*Pars::sigma/Pars::NU);
	printf("+-------------------------------------------------+\n");
	return;
}

void save_data::save_summary_log(Particle& par){
	// Cancel the saving procedure if flag is closed
	if (Pars::flag_save_parameter == false){return;}
	
	// Log saving before simulation
	std::ofstream data;
	data.open("output/Parameter.dat");
	data << "#=================================================#\n"
	     << "+---------------- SIMULATION LOG -----------------+\n"
		 << "#=================================================#\n\n";
	
	data << "+----------- Simulation Flow Parameter -----------+\n";
	data << std::fixed << std::setprecision(2)
	     << "Reynolds number (RE)             : "; data.width(8); data << std::right << Pars::RE    << " [-]" << "\n"
		 << "Freestream velocity (U)          : "; data.width(8); data << std::right << Pars::u_inf << " m/s"<< "\n"
 		 << "Fluid density (rho)              : "; data.width(8); data << std::right << Pars::RHO   << " kg/m^3" << "\n"
		 << "Fluid viscosity (nu)             : "; data.width(8); data << std::right << Pars::NU    << " m^2/s" << "\n"
		 << "+-------------------------------------------------+\n\n";

	data << "+----------- Simulation Setting Option -----------+\n"
	     << "Body option                             : " << "type " << Pars::opt_body << "\n"
		 << "Initialization option                   : " << "type " << Pars::opt_init << "\n"
		 << "Neighbor search option                  : " << "type " << Pars::opt_neighbor << "\n"
		 << "Penalization option                     : " << "type " << Pars::opt_pen << "\n";
	data << std::fixed << std::setprecision(4)
	     << "Number of save data                     : "; data.width(6); data << std::right << Pars::nt_data << "\n"
		 << "Saving step interval                    : "; data.width(6); data << std::right << Pars::step_inv << "\n"
		 << "Support radius factor                   : "; data.width(6); data << std::right << Pars::r_sup << "\n"
		 << "Buffer radius factor                    : "; data.width(6); data << std::right << Pars::r_buff << "\n"
		 << "Maximum resolution level                : "; data.width(6); data << std::right << Pars::max_level << "\n"
		 << "Core size                               : "; data.width(6); data << std::right << Pars::sigma << " m\n"
		 << "Time step                               : "; data.width(6); data << std::right << Pars::dt << " s\n";
	data << std::fixed << std::setprecision(2)
		 << "Total simulation time                   : "; data.width(6); data << std::right << Pars::sim_time << " s\n"
		 << "Total iteration step                    : "; data.width(6); data << std::right << Pars::nt << " \n"
		 << "+-------------------------------------------------+\n\n";

	data << "+----------- Simulation Parameter Data -----------+\n"
	     << "Domain x length                    :"; data.width(12); data << std::right << Pars::lxdom << " m\n"
		 << "Domain y length                    :"; data.width(12); data << std::right << Pars::lydom << " m\n"
		 << "Origin gap length                  :"; data.width(12); data << std::right << Pars::xdom  << " m\n"
		 << "Reference length (D)               :"; data.width(12); data << std::right << Pars::Df    << " m\n"
		 << "Plate_thick                        :"; data.width(12); data << std::right << Pars::Df*Pars::H_star   << " m\n";
	data << std::fixed << std::setprecision(4)
		 << "Courant number (C < 1.0)           :"; data.width(12); data << std::right << Pars::Courant << "\n"
		 << "Diffusion number (Phi < 0.5)       :"; data.width(12); data << std::right << Pars::Diffusion << "\n"
		 << "Stability Criteria (Re_h < 1.0)    :"; data.width(12); data << std::right << 100 * Pars::sigma*Pars::sigma/Pars::NU << "\n"
		 << "Turbulent scaling                  :"; data.width(12); data << std::right << std::ceil(1.0e0 / Pars::Courant) << "\n"
		 << "Maximum time step (Re_h criteria)  :"; data.width(12); data << std::right << 0.25*Pars::sigma*Pars::sigma/Pars::NU << " s\n"
		 << "+-------------------------------------------------+\n\n";
	data << "Number of initialized particle     :"; data.width(12);data << std::right << par.num << "\n";

	data << std::fixed << std::setprecision(-1);
	
	// End of writting simulation setting
	data.close();
	return;
}