#include "data_saving.hpp"

// Saving the body state
// BODY SURFACE
void save_data::save_state(const Body &b, string s){
	// Cancel the body saving if the flag turn off
	if (Pars::flag_save_body == false){return;}

	// Message Log
	std::cout << "\nSaving body data ...\n";

	// ! Saving parameters
	std::string name1;
	std::ofstream data;
	
	// Save the body node data
	name1.append("output/body_node_");
	name1.append(s);
	name1.append(".csv");
	data.open(name1.c_str());
	data << "node" << ","
		 << "xb" << "," 
		 << "yb" // << ","
		 // << "zb" // << ","
		 << "\n";

	for (int i = 0; i < b.num; i++)
	{
		data << i << ","
			 << b.x[i] << ","
			 << b.y[i] // << ","
			//  << b.z[i] << ","
			 << "\n";
	}
	data.close();

	// Save the body panel data
	name1 = "output/body_panel_";
	name1.append(s);
	name1.append(".csv");
	data.open(name1.c_str());
	data << "panel" << ","
		 << "x_mid" << "," 
		 << "y_mid" << "," 
		//  << "z_mid" << "," 
		 << "x_normal" << "," 
		 << "y_normal" // << "," 
		//  << "z_normal" << "," 
		//  << "nodeList" 
		 << "\n";

	for (int i = 0; i < b.n_panel; i++)
	{
		data << i << ","
			 << b.x_m[i] << "," 
			 << b.y_m[i] << "," 
			 //  << b.z_m[i] << "," 
			 << b.x_n[i] << "," 
			 << b.y_n[i] // << "," 
			 //  << b.z_n[i] << "," 
			 //  << "nodeList" //<< ","
			 << "\n";
	}
	data.close();
}

// PARTICLE DISTRIBUTION
void save_data::save_state(Particle &p, string s, int type)
{
	/* Particle distribution saving state description
       * @param p    : particle data to be saved
       * @param type : the type of data saving
         - Type 0 -> Save all particle data
		 - Type 1 -> Save the particle at the interested domain only
       * @param s    : the name for the state
    */
	
	// Store neighbor ID flag   || INPUT MANUALLY here (save the list of each par ID in verlet list)
	bool _ngh = false;
	bool _ngh_NUM = true;
	bool _pressure = Pars::flag_pressure_calc;

	// internal variables
	std::string name;
	std::ofstream ofs;
    int np = p.x.size();             // Number of particle
	double limdom = 1;          // Limitation for type 1 saving

	printf("\nSaving the particle state ...\n");
		
	// Save common data (information of free particles):
	name.append("output/particle_state_");
	name.append(s);
	name.append(".csv");
	ofs.open(name.c_str());
	ofs << "" << "xp" 
		<< "," << "yp" 
		<< "," << "gpz"
		<< "," << "vor" 
		<< "," << "up" 
		<< "," << "vp" 
		<< "," << "sp" 
		<< "," << "active" 
		<< "," << "chi" 
		// << "," << "R" 
		// << "," << "Basis_CELL_ID" 
		// << "," << "CELL_ID"
		;
	if (_pressure){
		ofs << "," << "pressure";
	}
	if (_ngh_NUM){
		ofs << "," << "ngh_num";
	}
	if (_ngh){
		ofs << "," << "ngh_ID";
	}
	ofs << "\n";

	for (int i = 0; i < np; i++)
	{
		if (type == 0)	// Save whole domain
		{
			ofs << "" << p.x[i]
				<< "," << p.y[i]
				<< "," << p.gz[i]
				<< "," << p.vorticity[i]
				<< "," << p.u[i] - Pars::ubody 
				<< "," << p.v[i] - Pars::vbody
				<< "," << p.s[i]
				<< "," << p.isActive[i]
				<< "," << p.chi[i]
				// << "," << p.R[i]
				// << "," << p.basis_label[i]
		    	// << "," << p.cell_label[i]
				;
			if (_pressure){
				ofs << "," << p.P[i];
			}
			if (_ngh_NUM){
				ofs << "," << p.neighbor[i].size();
			}
			if (_ngh){
				ofs << "," ;
				for (int j = 0; j < p.neighbor[i].size(); j++)
				{
					if(j == 0)
						ofs << p.neighbor[i][j];
					else
						ofs << "." << p.neighbor[i][j];
				}
			}
			ofs << "\n";
		}
		else if (type == 1)	// Only save the interest domain
		{
			if (p.x[i] > -limdom && p.x[i] < Pars::Df + limdom && p.y[i] > -limdom && p.y[i] < limdom)
			{
				ofs << "" << p.x[i]
					<< "," << p.y[i]
					<< "," << p.gz[i]
					<< "," << p.vorticity[i]
					<< "," << p.u[i] - Pars::ubody 
					<< "," << p.v[i] - Pars::vbody
					<< "," << p.s[i]
					<< "," << p.isActive[i]
					<< "," << p.chi[i]
					// << "," << p.R[i]
					// << "," << p.basis_label[i]
		            // << "," << p.cell_label[i]
					;
				if (_pressure){
					ofs << "," << p.P[i];
				}
				if (_ngh_NUM){
					ofs << "," << p.neighbor[i].size();
				}
				if (_ngh){
					ofs << "," ;
					for (int j = 0; j < p.neighbor[i].size(); j++)
					{
						if(j == 0)
							ofs << p.neighbor[i][j];
						else
							ofs << "." << p.neighbor[i][j];
					}
				}
				ofs << "\n";
			}
		}
	}
	ofs.close();
}
