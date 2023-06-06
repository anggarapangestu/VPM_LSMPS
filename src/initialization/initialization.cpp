#include "initialization.hpp"

initialization::initialization(/* args */)
{
}

initialization::~initialization()
{
}

// Particle Distribution Initialization Manager
void initialization::generate(const Body &b, Particle &p){
    printf("\n+------------- Particle Distribution -------------+\n");
    // ========== INITIALIZATION: Particles Generation ========== //
    // computational time accumulation
    clock_t t = clock();
    
    // generate the particle position distribution
    printf("Generating particle distribution ...\n");
	if(Pars::DIM == 2)
	{
		if (Pars::opt_init == 1) // Single Resolution Particle Distribution
        {
            printf("<+> Type 1: Single resolution\n");
            this->init_2d_single_res(p);
        }
        else if (Pars::opt_init == 2) // Two Level Resolution Single Block Particle Distribution
        {
            printf("<+> Type 2: Two level resolution single block\n");
            this->init_2d_multi_res_single_block(b, p);
        }
        // [Not Ready YET]
        else if (Pars::opt_init == 3) // Multi Resolution Multi Block Particle Distribution
        {
            printf("<+> Type 3: Multi resolution multi block\n");
            this->init_2d_multi_res_multi_block(b, p);
        }
        else if (Pars::opt_init == 4) // Multi Resolution Body Adjusted Particle Distribution
        {
            printf("<+> Type 4: Multi resolution body adjusted\n");
            this->init_2d_multi_res_body_adjusted(b, p);
        }
        else if (Pars::opt_init == 0) // [Test] Initialization
        {
            // initialization_step.init_particle(body, particle);	// <<<========= Type 1
		    // initialization_step.init_domain(particleBase);		// <<<========= Type 2
        }
	}
    else if(Pars::DIM == 3)
	{
		if (Pars::opt_init == 1) // Single Resolution Particle Distribution
        {
            printf("<+> Type 1: Single resolution\n");
            this->init_domain_3d(p);
        }
	}

    // particle generation initialization summary time display
    t = clock() - t;
    printf("<+> Number of initialize particle       : %8d \n", p.num);
	printf("<-> Particle initialization comp.\n");
    printf("    time:                              [%f s]\n", (double)t/CLOCKS_PER_SEC);

    
    // THE CODE BELOW MUST BE SEPARATED FROM THIS

    // ========== INITIALIZATION: Chi and Initial Active Particle Calculation ========== //
    // Adjust the computational time calculation
    t = clock();

    // Calculate the nearest particle distance to body surface
    d_geom.distance_calc(p,b);
    
    // Adjust the mollification parameter
    printf("Calculate chi parameter ...\n");
    p.chi.resize(p.num, 0.0e0);
    for (int i = 0; i < p.num; i++){
        // Active particle evaluation
        if (p.R[i] < Pars::init_act_dist) // Inside the active region, [adjusted MANUAL from global]
        {
            p.isActive[i] = true;
        }

        // Chi calculation and evaluation
        if (p.R[i] > Pars::hmollif) // Outside the body and mollification region
        {
            p.chi[i] = 0.0e0;
        }
        else if (std::abs(p.R[i]) <= Pars::hmollif) // Inside the mollification region (transition region)
        {
            p.chi[i] = -0.5e0 * (-1.0e0 + p.R[i] / Pars::hmollif + (1.0e0 / Pars::pi) * sin(Pars::pi * p.R[i] / Pars::hmollif));
        }
        else if (p.R[i] < Pars::hmollif) // Inside the body region only
        {
            p.chi[i] = 1.0e0;
        }
    }

    // chi calculation computational time display
    t = clock() - t;
	printf("<-> Chi calculation comp. time:        [%f s]\n", (double)t/CLOCKS_PER_SEC);
    

    printf("+-------------------------------------------------+\n");

}

// initialization_step.continue_simulation(body, particle, nt_start);   // Read the last particle data
void initialization::continue_simulation(const Body &b, Particle &par, int step){
    printf("\n+--------------- Particle Reading ----------------+\n");
    // ========== INITIALIZATION: Particles Reading ========== //
    clock_t t = clock();
    // Writting the file name from given step time
    printf("Reading particle data %d ...\n", step);
    std::string filename = "output/particle_state_";
    {
        // Digit of the maximum iteration
        int maxDigitLen = 1 + std::floor(std::log10(Pars::nt));
        
        // Digit of the current step
        int iterDigitLen = 1;
        if (step != 0){iterDigitLen = 1 + std::floor(std::log10(step));}
        
        // Procedure to get the name "particle_state_####.csv"
        int addDigit = maxDigitLen - iterDigitLen;
        for (int _spc = 0; _spc < addDigit; _spc++)
            filename.append("0");
        filename.append(std::to_string(step));
        filename.append(".csv");
    }
    
    // Read the data file
    std::fstream readData;
    readData.open(filename, std::ios::out | std::ios::in);

    // Internal variable
    char text[200];
    std::string value = "";
    double VAL;
    int first = 0;
    int loc;

    // Read all data line by line
    double xp,yp,gp,vor,up,vp,sp,chi;
    bool active;
    while(readData.peek()!=EOF){
        readData >> text;
        if (first == 0){first ++; continue;}
        loc = 0;

        // The sequence of the location index
        // xp,yp,gpz,vor,up,vp,sp,active
        // xp,yp,gpz,vor,up,vp,sp,active,chi,pressure,ngh_num
        // 0  1   2   3  4  5   6   7    8      9      10

        // Check for the line
        for (auto val:text){
            // Check the line terminator
            if (val == '\0'){
                VAL = std::stod(value);
                active = VAL;
                value = "";
                loc ++;
                break;
            }

            // End the value digit
            if (val == ','){
                VAL = std::stod(value);
                if (loc == 0) {xp = VAL;}
                else if (loc == 1) {yp = VAL;}
                else if (loc == 2) {gp = VAL;}
                else if (loc == 3) {vor = VAL;}
                else if (loc == 4) {up = VAL;}
                else if (loc == 5) {vp = VAL;}
                else if (loc == 6) {sp = VAL;}
                // else if (loc == 7) {active = VAL;}
                // else if (loc == 8) {chi = VAL;}
                value = "";
                loc ++;
                continue;
            }

            // Take all value of the string character
            value = value + val;
        }

        int _level = Pars::max_level - int(std::log2(sp/Pars::sigma));
        // Perform the geometry input data
        par.num ++;
        par.x.push_back(xp);
        par.y.push_back(yp);
        par.s.push_back(sp);
        par.level.push_back(_level);
        par.u.push_back(up);
        par.v.push_back(vp);
        par.vorticity.push_back(vor);
        par.gz.push_back(gp);
        // par.isActive.push_back(active);
        // par.chi.push_back(chi);
    }
    readData.close();

    // Resize other data
    par.isActive.clear();
    par.isActive.resize(par.num, false);
    
    // particle reading initialization summary time display
    t = clock() - t;
    printf("<+> Number of initialize particle       : %8d \n", par.num);
	printf("<-> Particle initialization comp.\n");
    printf("    time:                              [%f s]\n", (double)t/CLOCKS_PER_SEC);


    // ========== INITIALIZATION: Chi and Initial Active Particle Calculation ========== //
    // Adjust the computational time calculation
    t = clock();

    // Calculate the nearest particle distance to body surface
    d_geom.distance_calc(par,b);
    
    // Adjust the mollification parameter
    printf("Calculate chi parameter ...\n");
    par.chi.resize(par.num, 0.0e0);
    for (int i = 0; i < par.num; i++){
        // Active particle evaluation
        if (par.R[i] < Pars::init_act_dist) // Inside the active region, [adjusted MANUAL from global]
        {
            par.isActive[i] = true;
        }

        // Chi calculation and evaluation
        if (par.R[i] > Pars::hmollif) // Outside the body and mollification region
        {
            par.chi[i] = 0.0e0;
        }
        else if (std::abs(par.R[i]) <= Pars::hmollif) // Inside the mollification region (transition region)
        {
            par.chi[i] = -0.5e0 * (-1.0e0 + par.R[i] / Pars::hmollif + (1.0e0 / Pars::pi) * sin(Pars::pi * par.R[i] / Pars::hmollif));
        }
        else if (par.R[i] < Pars::hmollif) // Inside the body region only
        {
            par.chi[i] = 1.0e0;
        }
    }

    // chi calculation computational time display
    t = clock() - t;
	printf("<-> Chi calculation comp. time:        [%f s]\n", (double)t/CLOCKS_PER_SEC);
    printf("+-------------------------------------------------+\n");
}



/*
Note for Initialization:
1. <Must> add option for initialization type [DONE]
2. 
*/