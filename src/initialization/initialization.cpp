#include "initialization.hpp"
#include "../grid_block/generateGrid.hpp"
#include "../velocity_calculation/velocity_calc.hpp"
#include "../grid_block/gridNodeAdapt.hpp"
#include "../neighbor/neighbor.hpp"
#include "../adaptation/adaptation.hpp"

#include <bits/stdc++.h>

// ===================================================== //
// +------------------ Public Method ------------------+ //
// ===================================================== //
// #pragma region PUBLIC_METHOD

/**
 *  @brief  Particle initialization manager. Generate the particle based on the given
 *  criteria defined by user in global.cpp.
 *   GENERATE:   [1] Coordinate data (x,y,z); [2] Particle size (sigma); 
 *               [3] Tree hierarchy (level, nodeID) [An additional]
 *   CALCULATE:  [1] Flag Marker (body part, active, inside body, near surface);
 *               [2] Penalization (chi) {Upon body existance}
 *   INITIALIZE: [1] Velocity (u,v,w); [2] Vorticity (vortx, vorty, vortz);
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for particle generation.
 *  @param  _gridNode  Grid node data container used for particle generation.
 * 
 *  @return No return.
*/
void initialization::initialize_particle(Particle &_par, const std::vector<Body> &_BL, GridNode &_grid){
    /** PROCEDURE:
     *   [1] Generate the particle geometry data (initialized basic data)
     *        > Type 1: Start over (see generate particle subroutine)
     *        > Type 2: Read the particle data -> Generate grid node
     *   [2] Set up the body related flag (also calculate chi for penalization) [only if the body existed]
     *   [3] Update the active flag
    */
    
    printf("\n+------------ Particle Initialization ------------+\n");
    // Computational time accumulation
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // [1] Initialize the particle based on starting state option
    // **********************************************************
    // Start over initialization
    if (Pars::opt_start_state == 0){
        //  *generate the particle based on the body data
		this->generate_particle(_par, _BL, _grid);
	}
    
    // Resume the simulation from the given state
    else if(Pars::opt_start_state == 1){
        //  *read the particle data
        if (DIM == 2){
            this->read_2d_particle(_par, Pars::resume_step);
        }else if (DIM == 3){
            this->read_3d_particle(_par, Pars::resume_step);
        }

        // Generate grid node for type 5 initialization only
        if (Pars::opt_init_particle == 5){
            printf("%sGenerate grid block ... %s\n", FONT_CYAN, FONT_RESET);
            // Create grid from particle data list
            generateGrid grid_step;
            grid_step.createNode(_grid, _par);
        }
	}

    
    // [2] Evaluate all the body flag of the particle
    // **********************************************
    #if (N_BODY != 0)
        geometry geom_tool;
        geom_tool.eval_near_body(_par, _BL);    // Calculate the near body part ID
        geom_tool.eval_body_flag(_par, _BL);    // Calculate near surface flag, inside body flag, chi
    #endif


    // [3] Evaluate the active particle
    // ********************************
    // Initialize the active sign
    _par.isActive.resize(_par.num, false);

    if (Pars::opt_start_state == 0 && N_BODY != 0){
        // For start from begining
        #pragma omp parallel for
        for (int i = 0; i < _par.num; i++){
            // Only change the flag for near surface particle
            if (_par.isNearSurface[i] == true){
                _par.isActive[i] = true;
            }
        }
	}
    else if(Pars::opt_start_state == 1){
        // For resume simulation
        this->set_active_sign(_par, false);
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
    printf("<+> Number of initialize particles      : %8d \n", _par.num);
	printf("<-> Particle initialization computational\n");
    printf("    time:                              [%f s]\n", _time);

    printf("+-------------------------------------------------+\n");
    return;
}


/**
 *  @brief  Generate the initial value of vorticity field.
 *   CALCULATE: [1] Vorticity (a specific data); [2] Active Flag Marker, [3] Initial velocity
 *  @param  _particle  Particle data for vorticity calculation.
 *  @param  _remeshTool  The tool and member data for particle adaptation.
 *  @param  _baseGrid  The grid base for particle adaptation.
 *  @param  _bodyList  The list of body inside the domain.
 * 
 *  @return No return.
*/
void initialization::initialize_vorticity(Particle &_par, GridNode &_baseGrid, const std::vector<Body> &_bodyList){
    // Check the vorticity initialization type
    #if (N_BODY == 0)
        // Upon no body existance, the vorticity must be generated
        if (Pars::opt_init_vorticity == 0){
            ERROR_LOG << "Not enough simulation source! \nNo initialized vorticity neither initialized body!";
            throw std::runtime_error("[ERROR] An attempt to simulate a blank domain!");
        }
    #else
        std::cout << "Generate body adapted node ...\n";
        // Upon body existance, the vorticity may not be generated
        if (Pars::opt_init_vorticity == 0 || Pars::opt_start_state == 1){
            // Skip when initial vorticity is not needed
            // Also skip when a resume simulation is run
            return;
        }else{
            // Check whether need additional initialized vorticity
            WARNING_LOG << "An additional vorticity will be generated!\n";
            std::cout << "<!> Are you sure want to proceed (?) ("
			          << FONT_GREEN << "yes" << FONT_RESET << "/" 
			          << FONT_RED   << "no"  << FONT_RESET << ")\n" 
			          << "    Type here: ";
	
            // Prompt to continue run the simulation
            bool proceed = false; std::string _cmd; std::cin >> _cmd;
            if (_cmd == "yes" || _cmd == "Yes" || _cmd == "Y" || _cmd == "y") proceed = true;
            else proceed = false;
            if (!proceed) exit(1);
        }
    #endif

    // Display prompt
    printf("\n+----------- Vorticity Field Generation ----------+\n");

    // Computational time accumulation
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // [1] Vorticity field type selection
    // **********************************
    // Initial calculation
    this->calculate_vorticity(_par);
    this->set_active_sign(_par, false);
    
    // Particle refinement calculation (max 10 refinement)
    if (Pars::opt_init_particle == 5){
    for (int ctr = 0; ctr < 2*Pars::max_level; ctr++){
        // // Get the initial particle number
        // int i_num = _par.num;
        
        // Perform the resolution adaptation
        neighbor nghTools;
        nghTools.neigbor_search(_par, _baseGrid);

        Particle *tempPar = new Particle;
        *tempPar = _par;
        adaptation adapTools;
        bool _adapt = adapTools.get_adaptation(_par, tempPar, _baseGrid);

        if (_adapt){
            // Update the particle value
            _par = *tempPar;
            // Update the vorticity and active sign 
            this->calculate_vorticity(_par);
            this->set_active_sign(_par, false);
        }
        delete tempPar;

        // GridNodeAdapt adaptTools;
        // bool _adapt = adaptTools.get_adaptation(_baseGrid, _par, _par, -1);
        // if (_adapt){
        //     // Update the particle distribution
        //     Particle* _tempPar = new Particle;
        //     adaptTools.take_particle_pointer(_tempPar); // Take the new distribution here ...
        //     _par = *_tempPar;   // Update the particle here ...
        //     delete _tempPar;

        //     // Update the particle size (this is the drawback of current subroutine)
        //     for (int i = 0; i < _par.num; i++)
        //     _par.s[i] = Pars::sigma * Pars::intPow(2, Pars::max_level - _par.level[i]);
        //     // _par.s[i] = _baseGrid.gridSize/(Pars::intPow(2, _par.level[i]) * _baseGrid.baseParNum);

        //     // Update the vorticity and active sign 
        //     this->calculate_vorticity(_par);
        //     this->set_active_sign(_par, false);
        // }
        
        // // Get the initial particle number
        // int f_num = _par.num;
        
        // Check whether there are change in distribution number or not
        if (!_adapt/*f_num == i_num*/){
            std::cout << "Particle adaptation is well performed ...\n";
            break;
        }
    }    
    }

    // [2] Calculate the initial velocity
    // **********************************
    VelocityCalc velocity_tool;
    velocity_tool.get_velocity(_par, _baseGrid, 0);


    // Vorticity generation summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Vorticity generation comp. time    [%f s]\n", _time);

    printf("+-------------------------------------------------+\n");
    return;
}


/**
 *  @brief  Update the particle active sign by considering the value of particle
 *  vorticity.
 * 
 *  @param _particle The evaluated particle.
 *  @param _flag    The flag to reduce non-active particle.
 *  Set the vorticity property to 0 for non-active particle, when activated.
*/
void initialization::set_active_sign(Particle &par, bool flag){
    // **Find the maximum vorticity
    double vor_max = 0.0e0;
    for (int i = 0; i < par.num; i++){
        if (std::abs(par.vorticity[i]) > vor_max)
        vor_max = std::abs(par.vorticity[i]);
    }

    // **Evaluate the active flag particle
    par.isActive.clear(); par.isActive.resize(par.num,false);
    const double SIGNIFICANCE_LIMIT = Pars::active_sig*vor_max;
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        double vor = std::abs(par.vorticity[i]);
        par.isActive[i] = (vor >= SIGNIFICANCE_LIMIT) ? true : false;

        // Additional check (The particle insider near surface is categorized as active particle also)
        #if (N_BODY > 0)
            if (par.isNearSurface[i] == true){
                par.isActive[i] = true;
            }
        #endif 
    }
    

    // **Set the vorticity to zero for non active particle
    if (flag == true){
        #if (DIM == 2)
            // #pragma omp parallel for
            for (int i = 0; i < par.num; i++){
                if (par.isActive[i] == false){
                    // To prevent truncated error, set the value of vorticity to 0.0
                    par.vorticity[i] = 0.0e0;
                    par.gz[i] = 0.0e0;
                }
            }
            // // [OLD SECTION] not use but for reference
            // // Find the maximum absolute value of vorticity
            // double vor_max = 0.0e0;
            // for (int i = 0; i < par.num; i++){
            //     if (std::abs(par.vorticity[i]) > vor_max)
            //     vor_max = std::abs(par.vorticity[i]);
            // }

            // // Evaluate the active flag particle
            // par.isActive.clear();par.isActive.resize(par.num,false);
            // #pragma omp parallel for
            // for (int i = 0; i < par.num; i++){
            //     double vor = std::abs(par.vorticity[i]);
            //     if (vor >= Pars::active_sig*vor_max /* || par.isActive[i] == true*/){
            //         par.isActive[i] = true;
            //     }
            //     else{
            //         par.isActive[i] = false;
            //         // To prevent truncated error, set the value of vorticity to 0.0
            //         par.vorticity[i] = 0.0e0;
            //         par.gz[i] = 0.0e0;  
            //     }
            // }
        
        #elif (DIM == 3)
            // #pragma omp parallel for
            for (int i = 0; i < par.num; i++){
                if (par.isActive[i] == false){
                    // To prevent truncated error, set the value of vorticity to 0.0
                    par.vorticity[i] = 0.0e0;
                    par.vortx[i] = 0.0e0;
                    par.vorty[i] = 0.0e0;
                    par.vortz[i] = 0.0e0;
                }
            }
        #endif
    }
    return;
}


// #pragma endregion

// ===================================================== //
// +------------------ Private Method -----------------+ //
// ===================================================== //
// #pragma region PRIVATE_METHOD

/**
 *  @brief  Particle generator distribution manager. Generate the particle data
 *  using the selected type of distribution by user in global.cpp file.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for particle generation.
 *  @param  _gridNode  Grid node data container used for particle generation.
*/
void initialization::generate_particle(Particle &p, const std::vector<Body> &bL, GridNode &g){
    // Generate the particle position distribution
    printf("Generating particle data ...\n");

    /** The initialization section generates:
     *    [1] Coordinate data (x,y,z); 
     *    [2] Particle size (sigma) and level; 
     *    [3] Tree hierarchy (nodeID) [An additional]
    */

    #if (DIM == 2)  // Initialization for 2D simulation
        switch (Pars::opt_init_particle){
        case 0:
            // A Testing Initialization
            printf("%sType 0: Testing distribution %s\n", FONT_CYAN, FONT_RESET);

            //  * Put your initialization testing here ...
            // this->init_2d_test_1(p);        // Type 1 Initialization
            this->init_2d_test_2(p);        // Type 2 Initialization
            break;
        case 1:
            // Single Resolution Particle Distribution
            printf("%sType 1: Single resolution %s\n", FONT_CYAN, FONT_RESET);
            this->init_2d_single_res(p);
            break;
        case 2:
            // Two Level Resolution Single Block Particle Distribution
            printf("%sType 2: Two level resolution single block %s\n", FONT_CYAN, FONT_RESET);
            this->init_2d_multi_res_single_block(p, bL);
            break;
        case 3:     // [Not Ready YET]
            // Multi Resolution Multi Block Particle Distribution
            printf("%sType 3: Multi resolution multi block %s\n", FONT_CYAN, FONT_RESET);
            this->init_2d_multi_res_multi_block(p, bL);
            break;
        case 4:
            // Multi Resolution Body Adjusted Particle Distribution
            printf("%sType 4: Multi resolution body adjusted %s\n", FONT_CYAN, FONT_RESET);
            this->init_2d_multi_res_body_adjusted(p, bL);
            break;
        case 5:
            // Grid Block-Based Particle Initialization
            printf("%sType 5: Grid block particle generation %s\n", FONT_CYAN, FONT_RESET);
            this->init_par_grid_block(p, bL, g);
            break;
        default:
            break;
        }
    #elif(DIM == 3) // Initialization for 3D simulation
		switch (Pars::opt_init_particle){
        case 1:
            // Single Resolution Particle Distribution
            printf("%sType 1: Single resolution %s\n", FONT_CYAN, FONT_RESET);
            this->init_3d_single_res(p);
            break;
        case 2:
            // Two Level Resolution Single Block Particle Distribution
            printf("%sType 2: Two level resolution single block %s\n", FONT_CYAN, FONT_RESET);
            this->init_3d_multi_res_single_block(p, bL);
            break;
        case 4:
            // Multi Resolution Body Adjusted Particle Distribution
            printf("%sType 4: Multi resolution body adjusted %s\n", FONT_CYAN, FONT_RESET);
            this->init_3d_multi_res_body_adjusted(p, bL);
            break;
        case 5:
            // Grid Block-Based Particle Initialization
            printf("%sType 5: Grid block particle generation %s\n", FONT_CYAN, FONT_RESET);
            this->init_par_grid_block(p, bL, g);
            break;
        default:
            break;
        }
	#endif

    // Initialize the physical properties
    p.u.resize(p.num, Pars::u_inf);     // Velocity in x direction
    p.v.resize(p.num, Pars::v_inf);     // Velocity in y direction
    p.gz.resize(p.num, 0.0e0);          // Vortex strength in z direction
    p.vorticity.resize(p.num, 0.0e0);   // Vorticity scalar field
    
    if (Pars::flag_pressure_calc)
    p.P.resize(p.num, 0.0e0);           // Pressure field
    
    #if DIM == 3
    p.w.resize(p.num, Pars::w_inf);     // Velocity in z direction
    p.vortx.resize(p.num, 0.0e0);       // Vorticity strength in x direction
    p.vorty.resize(p.num, 0.0e0);       // Vorticity strength in y direction
    p.vortz.resize(p.num, 0.0e0);       // Vorticity strength in z direction
    #endif

    return;
}

/**
 *  @brief  Calculates the vorticity of particle distribution.
 *         
 *  @param  _particle  Particle data container.
*/
void initialization::calculate_vorticity(Particle &_par){
    // Calculate the vorticity of particle distribution
    printf("Calculating the vorticity ...\n");

    // [1] Vorticity field type selection
    // **********************************
    #if (DIM == 2)
        switch (Pars::opt_init_vorticity){
        case 1: // Eliptical vorticity
            printf("%sType 1: Eliptical vorticity ... %s\n", FONT_CYAN, FONT_RESET);
            this->eliptic_vorticity(_par);
            break;
        case 2: // Perlman vorticity
            printf("%sType 2: Perlman vorticity%s\n", FONT_CYAN, FONT_RESET);
            this->perlman_vorticity(_par);
            // this->perlman_velocity_solution(_par, 0);
            break;
        case 3: // Lamb Oseen vorticity
            printf("%sType 3: Lamb Oseen vorticity%s\n", FONT_CYAN, FONT_RESET);
            this->lamb_oseen_vorticity(_par);
            break;
        case 4: // Reserve ...
            printf("%sType 3: Reserved ... %s\n", FONT_CYAN, FONT_RESET);
            break;
        default:
            break;
        }
    #elif (DIM == 3)
        switch (Pars::opt_init_vorticity){
        case 1: // Ring vortex
            printf("%sType 1: Ring vortex ... %s\n", FONT_CYAN, FONT_RESET);
            this->ring_vortex(_par);
            break;
        case 2: // Reserve ...
            printf("%sType 2: Ring vortex oblique ... %s\n", FONT_CYAN, FONT_RESET);
            this->ring_vortex_oblique(_par);
            break;
        case 3: // Reserve ...
            printf("%sType 2: Reserved ... %s\n", FONT_CYAN, FONT_RESET);
            break;
        default:
            break;
        }
    #endif

    return;
}

/**
 *  @brief  Read the particle data in 2D simulation from system at the requested iteration time step.
 *  The read data are the basic properties: coordinate, size, and important physical properties.
 *  It will follows the pattern in data write method (see, save_data/state_data.cpp).
 *         
 *  @param  _particle   The particle storage to save the data value.
 *  @param  _iteration  Iteration step at the position data to be read.
 * 
 *  @headerfile initialization.hpp
*/
void initialization::read_2d_particle(Particle &_par, int _iter){
    // Print message log!!!
    printf("Reading particle data %d ...\n", _iter);

    // Reset the data
    _par.num = 0;
    _par.s.clear();
    _par.level.clear();
    _par.isActive.clear();

    _par.x.clear();
    _par.y.clear();
    
    _par.u.clear();
    _par.v.clear();

    _par.gz.clear();
    _par.vorticity.clear();
    
    if (Pars::flag_pressure_calc)
    _par.P.clear();
    
    // Construct the file directory name
    #if (DATA_INTERPOLATION == 0)
        std::string fileName = "input/resume_data/particle_state_";
    #elif (DATA_INTERPOLATION == 1)
        std::string fileName = "output/particle_state_";
    #endif
    simUtil util_tool;
    fileName += util_tool.saveName(_iter) + ".csv";
    
    // MESSAGE_LOG << "The filename " << fileName << " is about to read!\n";

    // Internal variable
    std::vector<double> dataList;
    std::string line;
    std::string entry;
    // Data storage
    double xp,yp;
    double up,vp;
    double vor,gz,P,sp;
    bool act;
    
    // Read the data file
    std::ifstream reader;
    reader.open(fileName);
    if(!reader.is_open()){
        ERROR_LOG << "The filename " << fileName << " is not existed!\n";
        throw std::runtime_error("Error: The intended file is not existed!");
    }
    getline(reader, line);  // Put out the first line (header) before reading
    
    // Read all data line by line
    while(reader.peek()!=EOF){
        getline(reader, line);      // Take the data line
        dataList.clear();           // Reset the list
        entry = "";                 // Reset the entry

        // The sequence of the location index for 2D
        // xp,yp,gpz,vor,up,vp,sp,active,chi,pressure,ngh_num,ngh_ID
        // 0  1   2   3  4  5   6   7    8      9      10      11

        // Get data per data in line, put into the data list
        for (const auto &dig : line){
            // Take the entry as the value
            if (dig == ','){
                double value = std::stod(entry);
                dataList.push_back(value);
                entry = "";
                continue;
            }

            // Add the digit into the entry
            entry += dig;
        }
        // Take the last data
        dataList.push_back(std::stod(entry));

        // Take out the data from list
        for (size_t i = 0; i < dataList.size(); i++){
            // Aliasing to the value
            const auto &value = dataList[i];
            int loc = i;
            
            // Put the data into particle
            switch (loc){
                case 0: xp  = value; break;
                case 1: yp  = value; break;
                case 2: gz  = value; break;
                case 3: vor = value; break;
                case 4: up  = value; break;
                case 5: vp  = value; break;
                case 6: sp  = value; break;
                case 7: act = value; break;
                case 9: P = value; break;
                default: break;
            }
        }

        // Calculate the particle level
        int lvl = Pars::max_level;
        int mul = (sp + Pars::sigma/100)/Pars::sigma;   // Particle size multiplier
        while(mul != 1){
            --lvl;
            mul/=2;
        }

        // Put the data into Perform the geometry input data
        _par.num++;
        _par.s.push_back(sp);
        _par.level.push_back(lvl);
        _par.isActive.push_back(act);
        
        _par.x.push_back(xp);
        _par.y.push_back(yp);
        
        _par.u.push_back(up);
        _par.v.push_back(vp);
        
        _par.gz.push_back(gz);
        _par.vorticity.push_back(vor);

        if (Pars::flag_pressure_calc)
        _par.P.push_back(P);

        // MESSAGE_LOG << "Done particle " << _par.num << " !\n";
        
    }
    reader.close();
    return;
}

/**
 *  @brief  Read the particle data in 3D simulation from system at the requested iteration time step.
 *  The read data are the basic properties: coordinate, size, and important physical properties.
 *  It will follows the pattern in data write method (see, save_data/state_data.cpp).
 *         
 *  @param  _particle   The particle storage to save the data value.
 *  @param  _iteration  Iteration step at the position data to be read.
 * 
 *  @headerfile initialization.hpp
*/
void initialization::read_3d_particle(Particle &_par, int _iter){
    // Print message log!!!
    printf("Reading particle data %d ...\n", _iter);

    // Reset the data
    _par.num = 0;
    _par.s.clear();
    _par.level.clear();
    _par.isActive.clear();
    
    _par.x.clear();
    _par.y.clear();
    _par.z.clear();

    _par.u.clear();
    _par.v.clear();
    _par.w.clear();

    _par.vortx.clear();
    _par.vorty.clear();
    _par.vortz.clear();
    _par.vorticity.clear();
    
    if (Pars::flag_pressure_calc)
    _par.P.clear();
    
    // Construct the file directory name
    #if (DATA_INTERPOLATION == 0)
    std::string fileName = "input/resume_data/particle_state_";
    #elif (DATA_INTERPOLATION == 1)
    std::string fileName = "output/particle_state_";
    #endif
    simUtil util_tool;
    fileName += util_tool.saveName(_iter) + ".csv";

    // Internal variable
    std::vector<double> dataList;
    std::string line;
    std::string entry;
    // Data storage
    double xp,yp,zp;
    double up,vp,wp;
    double vortx,vorty,vortz;
    double vor,P,sp;
    bool act;

    // Read the data file
    std::ifstream reader;
    reader.open(fileName);
    if(!reader.is_open()){
        ERROR_LOG << "The filename " << fileName << " is not existed!\n";
        throw std::runtime_error("Error: The intended file is not existed!");
    }
    getline(reader, line);  // Put out the first line (header) before reading
    
    // Read all data line by line
    while(reader.peek()!=EOF){
        getline(reader, line);      // Take the data line
        dataList.clear();           // Reset the list
        entry = "";                 // Reset the entry

        // The sequence of the location index for 3D
        // xp,yp,zp,vortx,vorty,vortz,vor,up,vp,wp,sp,active,chi,pressure,ngh_num,ngh_ID
        // 0  1  2    3     4     5    6  7  8  9  10   11   12     13       14     15

        // Get data per data in line, put into the data list
        for (const auto &dig : line){
            // Take the entry as the value
            if (dig == ','){
                double value = std::stod(entry);
                dataList.push_back(value);
                entry = "";
                continue;
            }

            // Add the digit into the entry
            entry += dig;
        }
        // Take the last data
        dataList.push_back(std::stod(entry));

        // Take out the data from list
        for (size_t i = 0; i < dataList.size(); i++){
            // Aliasing to the value
            const auto &value = dataList[i];
            int loc = i;
            
            // Put the data into particle
            switch (loc){
                case 0: xp = value; break;
                case 1: yp = value; break;
                case 2: zp = value; break;
                case 3: vortx = value; break;
                case 4: vorty = value; break;
                case 5: vortz = value; break;
                case 6: vor = value; break;
                case 7: up = value; break;
                case 8: vp = value; break;
                case 9: wp = value; break;  
                case 10: sp = value; break;
                case 11: act = value; break;
                case 13: P = value; break;
                default: break;
            }
        }

        // Calculate the particle level
        int lvl = Pars::max_level;
        int mul = (sp + Pars::sigma/100)/Pars::sigma;   // Particle size multiplier
        while(mul != 1){
            --lvl;
            mul/=2;
        }

        // Put the data into Perform the geometry input data
        _par.num++;
        _par.s.push_back(sp);
        _par.level.push_back(lvl);
        _par.isActive.push_back(act);

        _par.x.push_back(xp);
        _par.y.push_back(yp);
        _par.z.push_back(zp);

        _par.u.push_back(up);
        _par.v.push_back(vp);
        _par.w.push_back(wp);

        _par.vortx.push_back(vortx);
        _par.vorty.push_back(vorty);
        _par.vortz.push_back(vortz);
        _par.vorticity.push_back(vor);

        if (Pars::flag_pressure_calc)
        _par.P.push_back(P);
        
    }
    reader.close();
    
    return;
}

// #pragma endregion

// ===================================================== //
// +---------------- Grid Block Method ----------------+ //
// ===================================================== //
// #pragma region GRID_BLOCK_METHOD

/** Exactly cut the particle outside the defined domain
 *   [0] : Leave as it is
 *   [1] : Trim the particle outside the defined domain
*/
#define TRIM_DOMAIN_FLAG 0

/** The flag to add boundary particle
 *   [0] : No boundary particle
 *   [1] : Add boundary particle
*/

#define ADD_BOUNDARY_PARTICLE_FLAG 0
/**
 *  @brief  Particle generator using grid block method.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for particle generation.
 *  @param  _gridNode  Grid node data container used for particle generation.
*/
void initialization::init_par_grid_block(Particle &par, const std::vector<Body> &_bodyList, GridNode &_gridNode){
    // Generate Node List
    generateGrid grid_step;
    if (_bodyList.size() == 0){
        // Generating a node with single resolution
        std::cout << "Generate max level resolution node ...\n";
        grid_step.nodeGeneration(_gridNode);
    }else{
        // Generating a node body adapted
        std::cout << "Generate body adapted node ...\n";
        grid_step.nodeGeneration(_gridNode, _bodyList);
    }

    // Reserve the memory for particle data
    par.x.clear();
    par.y.clear();
    par.z.clear();
    par.nodeID.clear();
    par.level.clear();
    par.s.clear();

    // Initialize the divisor
    int parID = 0;          // Temporary particle ID local to node
    int _parNum = Pars::intPow(_gridNode.baseParNum, DIM);  // Number of particle inside node
    int _div[DIM];          // Particle ID divisor
    int _parIndex[DIM];     // Temporary index position of particle local to node
    
    // Initialize the divisor
    _div[0] = 1;
    for (int i = 1; i < DIM; i++){
        _div[i] = _div[i-1] * _gridNode.baseParNum;
    }

    // Generate particle by loop through all nodes
    for (auto &[_ID, _node] : _gridNode.nodeMap){
        // Only do particle generation to leaf node
        if (_node->isLeaf){
            #if TRIM_DOMAIN_FLAG == 1
                // Calculate shared properties of particle inside the current node
                double _parSize = _node->length / _gridNode.baseParNum;
                if (_node->isBoundary == true){
                    // Generate all particle inside the curren Node
                    for (int _locID = 0; _locID < _parNum; _locID++){
                        // Calculate local index coordinate inside the node from local ID
                        basis_loop(d) _parIndex[d] = (_locID/_div[d]) % _gridNode.baseParNum;

                        // Calculate all particle position
                            double _x = _node->pivCoor[0] + (0.5 + _parIndex[0])*_parSize;
                            if (_x > _gridNode.maxDomBound[0] || _x < _gridNode.minDomBound[0]) continue;
                        #if (DIM > 1)
                            double _y = _node->pivCoor[1] + (0.5 + _parIndex[1])*_parSize;
                            if (_y > _gridNode.maxDomBound[1] || _y < _gridNode.minDomBound[1]) continue;
                        #endif
                        #if (DIM > 2)
                            double _z = _node->pivCoor[2] + (0.5 + _parIndex[2])*_parSize;
                            if (_z > _gridNode.maxDomBound[2] || _z < _gridNode.minDomBound[2]) continue;
                        #endif
                        
                        // Update the location
                            par.x.push_back(_x);
                        #if (DIM > 1)
                            par.y.push_back(_y);
                        #endif
                        #if (DIM > 2)
                            par.z.push_back(_z);
                        #endif

                        // Assign other data
                        par.nodeID.push_back(_node->nodeID);
                        par.level.push_back(_node->level);
                        par.s.push_back(_parSize);

                        // Insert the current particle ID into the corresponding node
                        _node->parList.push_back(parID++);
                    }
                }else{
                    // Generate all particle inside the curren Node
                    for (int _locID = 0; _locID < _parNum; _locID++){
                        // Calculate local index coordinate inside the node from local ID
                        basis_loop(d) _parIndex[d] = (_locID/_div[d]) % _gridNode.baseParNum;

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
            #else
                // Calculate shared properties of particle inside the current node
                double _parSize = _node->length / _gridNode.baseParNum;

                // Generate all particle inside the curren Node
                for (int _locID = 0; _locID < _parNum; _locID++){
                    // Calculate local index coordinate inside the node from local ID
                    basis_loop(d) _parIndex[d] = (_locID/_div[d]) % _gridNode.baseParNum;

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
            #endif
        }
    }

    // Update the particle number
    par.num = parID;
    return;
}

// #pragma endregion

// ===================================================== //
// +----------------- Laplace Package -----------------+ //
// ===================================================== //
// #pragma region ADDITIONAL_LAPLACE_PACKAGE

// Boundary value
#define BC_LEFT     1       // Left type boundary condition [1:=dirichlet, 2:=neumann]
#define BC_RIGHT    1       // Right type boundary condition [1:=dirichlet, 2:=neumann]
#define BC_TOP      1       // Top type boundary condition [1:=dirichlet, 2:=neumann]
#define BC_BOTTOM   1       // Bottom type boundary condition [1:=dirichlet, 2:=neumann]
#define VOR_SIGNIFICANCE 1e-4  // The significant value of vorticity

#define GENERATE_RANDOM_SHIFT 1  // The flag to generate a random shift into the calculation

#define VELOCITY_BASED_POISSON 0  // The flag to generate a random shift into the calculation

/**
 *  @brief  [Dummy function]
*/
void initialization::initialize_laplace(Particle &_par){
    /** Procedure:
     *   [1] Initialize the boundary condition with label below
     *        -2:= negative y boundary (bottom)
     *         2:= positive y boundary (top)
     *        -1:= negative x boundary (left)
     *         1:= positive x boundary (right)  
     *   [2] Generate the particle inside domain
    */

    // Laplace parameter 
    double V_0 = 10.0;       // The potential at the left side (x=-sx)
    
    // PROCEDURE 1!
    // ***********
    // Set the location of each boundary
    // Illustration
    //    ___
    //  |     | 
    //  | ___ |
    // Each boundary is not necessary close the corner

    // Internal variable
    double minLoc[DIM];     // The minimum location of x and y
    double maxLoc[DIM];     // The maximum location of x and y
    double pivot[DIM];      // The pivot of each boundary

    // Assign the value
    minLoc[0] = -0.5 * Pars::lxdom;
    maxLoc[0] = 0.5 * Pars::lxdom;
    minLoc[1] = -0.5 * Pars::lydom;
    maxLoc[1] = 0.5 * Pars::lydom;

    // Generate the boundary value
    int bndCntX = std::ceil(Pars::lxdom / Pars::sigma)+1;    // The number (count) of boundary particle at x side
    int bndCntY = std::ceil(Pars::lydom / Pars::sigma)-1;    // The number (count) of boundary particle at y side
    double spaceX = Pars::lxdom / (bndCntX-1);
    double spaceY = Pars::lydom / (bndCntY+1);
    // pivot[0] = -0.5*(bndCntX-1)*spaceX;
    // pivot[1] = -0.5*(bndCntY+1)*spaceY;
    pivot[0] = minLoc[0];
    pivot[1] = minLoc[1];

    // Generate the boundary particle
    // Bottom part (labeled as -2)
    _par.num += bndCntX;
    _par.y.resize(_par.num, minLoc[1]);
    _par.boundaryLoc.resize(_par.num, -2);
    _par.boundaryVal.resize(_par.num, 0.0);
    for (int i = 0; i < bndCntX; i++){
        _par.x.push_back(pivot[0] + i*spaceX);
        // _par.y.push_back(minLoc[1]);
        // _par.boundaryLoc.push_back(-2);
        // _par.boundaryVal.push_back(BC_BOTTOM);
    }

    // Top part (labeled as 2)
    _par.num += bndCntX;
    _par.y.resize(_par.num, maxLoc[1]);
    _par.boundaryLoc.resize(_par.num, 2);
    _par.boundaryVal.resize(_par.num, 0.0);
    for (int i = 0; i < bndCntX; i++){
        _par.x.push_back(pivot[0] + i*spaceX);
        // _par.y.push_back(maxLoc[1]);
        // _par.boundaryLoc.push_back(2);
        // _par.boundaryVal.push_back(BC_TOP);
    }
    
    // Left part (labeled as -1)
    _par.num += bndCntY;
    _par.x.resize(_par.num, minLoc[0]);
    _par.boundaryLoc.resize(_par.num, -1);
    _par.boundaryVal.resize(_par.num, 0);
    for (int i = 0; i < bndCntY; i++){
        // _par.x.push_back(minLoc[0]);
        _par.y.push_back(pivot[1] + (1+i)*spaceY);
        // _par.boundaryLoc.push_back(-1);
        // _par.boundaryVal.push_back(BC_LEFT);
    }

    // Right part (labeled as 1)
    _par.num += bndCntY;
    _par.x.resize(_par.num, maxLoc[0]);
    _par.boundaryLoc.resize(_par.num, 1);
    _par.boundaryVal.resize(_par.num, V_0);
    for (int i = 0; i < bndCntY; i++){
        // _par.x.push_back(maxLoc[0]);
        _par.y.push_back(pivot[1] + (1+i)*spaceY);
        // _par.boundaryLoc.push_back(1);
        // _par.boundaryVal.push_back(BC_RIGHT);
    }

    
    // PROCEDURE 2!
    // ***********
    // Calculate the number of particle in each basis direction
    bndCntX = std::ceil(Pars::lxdom / Pars::sigma);    // The number (count) of boundary particle at x side
    bndCntY = std::ceil(Pars::lydom / Pars::sigma);    // The number (count) of boundary particle at y side
    pivot[0] = -0.5*bndCntX*Pars::sigma;
    pivot[1] = -0.5*bndCntY*Pars::sigma;
    int nTot = bndCntX * bndCntY;

    // Update the other properties
    _par.num += nTot;
    _par.s.resize(_par.num, Pars::sigma);
    _par.boundaryLoc.resize(_par.num, 0);
    _par.boundaryVal.resize(_par.num, 0.0);

    // Initialization for randomization shifting
    srand(time(0));
    double Rs = 0.5;

    // Generate the inner side particle
    for(int j = 0; j < bndCntY; j++){
        for(int i = 0; i < bndCntX; i++){
            double x = pivot[0] + (0.5+i)*Pars::sigma;
            double y = pivot[1] + (0.5+j)*Pars::sigma;

            #if (GENERATE_RANDOM_SHIFT == 1)
                double phi = ((double)std::rand()) / RAND_MAX;
                double mul = ((double)std::rand()) / RAND_MAX;
                double addX = std::cos(2*M_PI*phi) * mul * Rs * Pars::sigma;
                double addY = std::sin(2*M_PI*phi) * mul * Rs * Pars::sigma;
                x += addX;
                y += addY;
                if (x < minLoc[0] || x > maxLoc[0]){
                    x -= addX;
                }
                if (y < minLoc[1] || y > maxLoc[1]){
                    y -= addY;
                }
            #endif

            _par.x.push_back(x);
            _par.y.push_back(y);
        }
    }
    return;
}

/**
 *  @brief  Update the flag of particle near the domain boundary.
 *  NOTE: Actually generating the boundary condition.
 *         
 *  @param  _particle  Particle data container.
*/
void initialization::update_domain_boundary(Particle &_par){
    // Prompt display
    std::cout << "Updating the domain boundary ...\n";

    // Flag to calculate the velocity based on velocity poisson not the stream function poisson
    bool velPoisFlag = true;

    // Reserve the flag container data
    // _par.isBoundary = std::vector<bool>(_par.num, false);    // Update the flag size
    _par.boundaryLoc = std::vector<int>(_par.num, 0);           // Update the boundary location
    #if (VELOCITY_BASED_POISSON == 0)
        _par.boundaryVal = std::vector<double>(_par.num, 0.0);  // Update the boundary value
    #elif (VELOCITY_BASED_POISSON == 1)
        _par.vortx = std::vector<double>(_par.num, 0.0);  // Update the boundary value (FOR VELOCITY BC)
        _par.vorty = std::vector<double>(_par.num, 0.0);  // Update the boundary value (FOR VELOCITY BC)
    #endif

    // Get the domain extremes
    std::vector<double> maxPos(DIM, 0.0);
    std::vector<double> minPos(DIM, 0.0);

    for (int i = 0; i < _par.num; i++){
        // Update the x basis
        if (_par.x[i] > maxPos[0]) maxPos[0] = _par.x[i];
        if (_par.x[i] < minPos[0]) minPos[0] = _par.x[i];

        // Update the y basis
        if (_par.y[i] > maxPos[1]) maxPos[1] = _par.y[i];
        if (_par.y[i] < minPos[1]) minPos[1] = _par.y[i];

        #if (DIM == 3)
            // Update the z basis
            if (_par.y[i] > maxPos[1]) maxPos[1] = _par.y[i];
            if (_par.y[i] < minPos[1]) minPos[1] = _par.y[i];
        #endif
    }

    // Set the boundary value for each particle
    for (int i = 0; i < _par.num; i++){
        // Check the location of particle
        if (_par.x[i] < (minPos[0] + _par.s[i]/2.0)){
            // LEFT boundary
            _par.boundaryLoc[i] = -1;

            #if (BC_LEFT == 1)
                // A Diriclet boundary condition
                #if (VELOCITY_BASED_POISSON == 0)
                    _par.boundaryVal[i] = _par.phi_a[i];
                #elif (VELOCITY_BASED_POISSON == 1)
                    _par.vortx[i] = _par.u_a[i];
                    _par.vorty[i] = _par.v_a[i];
                #endif
            #elif (BC_LEFT == 2)
                // A Neumann boundary condition
                #if (VELOCITY_BASED_POISSON == 0)
                    _par.boundaryVal[i] = -_par.gy[i];
                #elif (VELOCITY_BASED_POISSON == 1)
                    _par.vortx[i] = _par.dudx[i];
                    _par.vorty[i] = _par.dvdx[i];
                #endif
            #endif
        }
        else if(_par.x[i] > (maxPos[0] - _par.s[i]/2.0)){
            // RIGHT boundary
            // if (std::abs(_par.vorticity[i]) < VOR_SIGNIFICANCE){
                _par.boundaryLoc[i] = 1;
                
                #if (BC_RIGHT == 1)
                    // A Diriclet boundary condition
                    #if (VELOCITY_BASED_POISSON == 0)
                        _par.boundaryVal[i] = _par.phi_a[i];
                    #elif (VELOCITY_BASED_POISSON == 1)
                        _par.vortx[i] = _par.u_a[i];
                        _par.vorty[i] = _par.v_a[i];
                    #endif
                #elif (BC_RIGHT == 2)
                    // A Neumann boundary condition
                    #if (VELOCITY_BASED_POISSON == 0)
                        _par.boundaryVal[i] = -_par.gy[i];
                    #elif (VELOCITY_BASED_POISSON == 1)
                        _par.vortx[i] = _par.dudx[i];
                        _par.vorty[i] = _par.dvdx[i];
                    #endif
                #endif
            // }
        }
        else if (_par.y[i] > (maxPos[1] - _par.s[i]/2.0)){
            // TOP boundary
            // if (std::abs(_par.vorticity[i]) < VOR_SIGNIFICANCE){
                _par.boundaryLoc[i] = 2;
                
                #if (BC_TOP == 1)
                    // A Diriclet boundary condition
                    #if (VELOCITY_BASED_POISSON == 0)
                        _par.boundaryVal[i] = _par.phi_a[i];
                    #elif (VELOCITY_BASED_POISSON == 1)
                        _par.vortx[i] = _par.u_a[i];
                        _par.vorty[i] = _par.v_a[i];
                    #endif
                #elif (BC_TOP == 2)
                    // A Neumann boundary condition
                    #if (VELOCITY_BASED_POISSON == 0)
                        _par.boundaryVal[i] = _par.gx[i];
                    #elif (VELOCITY_BASED_POISSON == 1)
                        _par.vortx[i] = _par.dudy[i];
                        _par.vorty[i] = _par.dvdy[i];
                    #endif
                #endif
            // }
        }
        else if( _par.y[i] < (minPos[1] + _par.s[i]/2.0)){
            // BOTTOM boundary
            // if (std::abs(_par.vorticity[i]) < VOR_SIGNIFICANCE){
                _par.boundaryLoc[i] = -2;
                
                #if (BC_BOTTOM == 1)
                    // A Diriclet boundary condition
                    #if (VELOCITY_BASED_POISSON == 0)
                        _par.boundaryVal[i] = _par.phi_a[i];
                    #elif (VELOCITY_BASED_POISSON == 1)
                        _par.vortx[i] = _par.u_a[i];
                        _par.vorty[i] = _par.v_a[i];
                    #endif
                #elif (BC_BOTTOM == 2)
                    // A Neumann boundary condition
                    #if (VELOCITY_BASED_POISSON == 0)
                        _par.boundaryVal[i] = _par.gx[i];
                    #elif (VELOCITY_BASED_POISSON == 1)
                        _par.vortx[i] = _par.dudy[i];
                        _par.vorty[i] = _par.dvdy[i];
                    #endif
                #endif
            // }
        }
    }

    std::cout << "Done updating domain boundary ...\n";

    return;
}

// #pragma endregion