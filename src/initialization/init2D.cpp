#include "initialization.hpp"

// =========================================
// ========= Particle Distribution =========
// =========================================

// [DONE] Uniform single resolution distribution
void initialization::init_2d_single_res(Particle &par)
{   
    // internal temporary variable
    double _y, _x;
    int i, j, count;

    // number of particle in x and y direction
    int _nx = std::ceil(Pars::lxdom / Pars::sigma);     // Number of particle in x direction
    int _ny = std::ceil(Pars::lydom / Pars::sigma);     // Number of particle in y direction
    int _nparticle = _nx * _ny;                 // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);

    // Generate the particle position
    count = 0;
    for (i = 0; i < _nx; i++)
    {
        for (j = 0; j < _ny; j++)
        {
            _x = -Pars::xdom + (i + 0.5)*Pars::sigma;
            _y = -Pars::sigma * _ny * 0.5 + (j + 0.5)*Pars::sigma;
            par.x[count] = _x;
            par.y[count] = _y;
            count ++;
        }
    }
    
    // Assign other particle properties
    par.num = _nparticle;
    par.s.resize(_nparticle, Pars::sigma);
    par.level.resize(_nparticle, Pars::max_level);
}

// [DONE] Two level resolution with finer resolution particle in single block distribution
void initialization::init_2d_multi_res_single_block(const Body &b, Particle &par)
{
    // internal temporary variable
    double _y, _x;
    double *bd_pos_max = new double[2]; // The maximum position for boundary x,y
    double *bd_pos_min = new double[2]; // The minimum position for boundary x,y
    int i, j, count;

    // ====== PARTICLE DISCRETIZATION GENERATION ====== //
    // number of particle in x and y direction
    int _nx = std::ceil(0.5 * Pars::lxdom / Pars::sigma);     // Number of particle in x direction
    int _ny = std::ceil(0.5 * Pars::lydom / Pars::sigma);     // Number of particle in y direction
    int _nparticle = _nx * _ny;                 // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);
    par.s.resize(_nparticle, 2 * Pars::sigma);
    par.level.resize(_nparticle, Pars::max_level - 1);
    par.num = _nparticle;

    // Generate the particle position
    count = 0;
    for (i = 0; i < _nx; i++)
    {
        for (j = 0; j < _ny; j++)
        {
            _x = -Pars::xdom + (2*i + 1.0)*Pars::sigma;
            _y = -Pars::sigma*_ny + (2*j + 1.0)*Pars::sigma;
            par.x[count] = _x;
            par.y[count] = _y;
            count ++;
        }
    }
    
    // ====== DIVIDE AND CONQUEROR PROCEDURE ====== //
    // define the boundary of finer resolution block
    bd_pos_max[0] = b.x[0]; // maximum x position of body node
    bd_pos_max[1] = b.y[0]; // maximum y position of body node
    bd_pos_min[0] = b.x[0]; // minimum x position of body node
    bd_pos_min[1] = b.y[0]; // minimum y position of body node
    
    // define the position limit of body surface node
    for (i = 1; i < b.num; i++){
        bd_pos_max[0] = bd_pos_max[0] > b.x[i] ? bd_pos_max[0] : b.x[i];
        bd_pos_max[1] = bd_pos_max[1] > b.y[i] ? bd_pos_max[1] : b.y[i];
        bd_pos_min[0] = bd_pos_min[0] < b.x[i] ? bd_pos_min[0] : b.x[i];
        bd_pos_min[1] = bd_pos_min[1] < b.y[i] ? bd_pos_min[1] : b.y[i];
    }

    // Put extension from the body limit to the finer resolution block limit (MANUAL EDIT)
    bd_pos_min[0] -= 0.75e0; // min x: Upstream
    bd_pos_max[0] += 1.50e0; // max x: Downstream
    bd_pos_min[1] -= 0.75e0; // min y: Lower body
    bd_pos_max[1] += 0.75e0; // max y: Upper body

    // Perform the devide and conqueror if the point inside the finer resolution block
    int* DnC_x = new int[4]{1,-1,-1,1};   // Divide and conqueror sum for x
    int* DnC_y = new int[4]{1,1,-1,-1};   // Divide and conqueror sum for y
    
    for (i = 0; i < par.num; i++){
        if ((par.x[i] < bd_pos_max[0] && par.x[i] > bd_pos_min[0]) && 
            (par.y[i] < bd_pos_max[1] && par.y[i] > bd_pos_min[1]))
        // DIVIDE AND CONQUEROR
        {
            // The particle is inside the finer resolution block
            _x = par.x[i];
            _y = par.y[i];
            
            // Point 1
            par.x[i] = _x + 0.5 * Pars::sigma;
            par.y[i] = _y + 0.5 * Pars::sigma;
            par.s[i] = Pars::sigma;
            par.level[i] = Pars::max_level;

            // Point 2 - 4
            for (j = 1; j < 4; j++){
                par.x.emplace_back(_x + 0.5 * Pars::sigma * DnC_x[j]);
                par.y.emplace_back(_y + 0.5 * Pars::sigma * DnC_y[j]);
                par.s.emplace_back(Pars::sigma);
                par.level.emplace_back(Pars::max_level);
                _nparticle ++;
            }
        }
    }
    
    // Deallocate the data
    delete DnC_x, DnC_y;
    delete bd_pos_max, bd_pos_min;

    // Assign other particle properties
    par.num = _nparticle;
}

// Multiresolution based on multiblock distribution
void initialization::init_2d_multi_res_multi_block(const Body &b, Particle &par){
    // TO BE CONTINUE
}
/*
void initialization::init_2d_multi_res_multi_block(const Body &b, Particle &par){
    // Calculate L_max, level_max, level_min
    double L_max;       // The factor size from the object size
    double max_grid_size;       // The factor size from the object size
    
    L_max = Pars::grid_factor * Pars::lx;
    this->level_max = 1 + std::round(std::log10(L_max / Pars::sigma)/std::log10(2));
    this->level_min = this->level_max + 1 - Pars::levels;

    max_grid_size = Pars::sigma * std::pow(2, this->level_max - 1);

    // Generate the base grid, set the level to min_level
    // Find the boundary position of the particle
    xmin_par[0] = *std::min_element(p_transformed.x.begin(), p_transformed.x.end()) + max_grid_size/2;
    xmin_par[1] = *std::min_element(p_transformed.y.begin(), p_transformed.y.end()) + max_grid_size/2;
    xmax_par[0] = *std::max_element(p_transformed.x.begin(), p_transformed.x.end());
    xmax_par[1] = *std::max_element(p_transformed.y.begin(), p_transformed.y.end());


    // Calculte the number of grid base particle in x and y direction
    for (size_t i = 0; i < 2; i++)
    {
        npar[i] = std::floor(std::abs(xmax_par[i] - xmin_par[i]) / max_grid_size);
    }

    std::vector<double> _xgrid;
    std::vector<double> _ygrid;
    
    // Create the grid base particle position
    for (size_t i = 0; i < npar[0]; i++)    // x direction
    {
        double _xBoxPos = xmin_par[0] + static_cast<double>(i) * max_grid_size;
        _xgrid.push_back(_xBoxPos);
    }
    
    for (size_t i = 0; i < npar[1]; i++)    // y direction
    {
        double _yBoxPos = xmin_par[1] + static_cast<double>(i) * max_grid_size;
        _ygrid.push_back(_yBoxPos);
    }

    for (size_t i = 0; i < npar[0]; i++)
    {
        for (size_t j = 0; j < npar[1]; j++)
        {
            this->GridBase.x.push_back(_xgrid[i]);
            this->GridBase.y.push_back(_ygrid[j]);
            this->GridBase.s.push_back(max_grid_size);
            this->GridBase.gz.push_back(0);
            //this->GridBase.vorticity.push_back(0);
            this->GridBase.level.push_back(this->level_min);
            
        }
    }
    
    this->GridBase.num = this->GridBase.x.size();   // Update the GridBase size
    _num = this->GridBase.num;
    
    // Find the grid neighbor
    std::vector<std::vector<int>> GridNeighbor = d_base_remeshing.inter_search_base_grid(this->GridBase);         // neighbor of base grid
    this->GridBase.neighbor = GridNeighbor;
    
    // ========================================
    // Calculate the body panel
    // Calculate the body boundary

    std::vector<double> _minBodyCoordinate(2,0); // center coordinate of obstacle
    std::vector<double> _maxBodyCoordinate(2,0); // center coordinate of obstacle

    std::vector<double> _centerPanelX;
    std::vector<double> _centerPanelY;
    std::vector<double> _normalPanelX;
    std::vector<double> _normalPanelY;

    double _nPanel = b.x.size() - 1;
    _centerPanelX.resize(_nPanel);
    _centerPanelY.resize(_nPanel);
    _normalPanelX.resize(_nPanel);
    _normalPanelY.resize(_nPanel);

    // Transform into circular
    Particle b_transformed; // Base grid particle (final interpolated coordinate)
    b_transformed.x.resize(b.x.size(), 0.0);
    b_transformed.y.resize(b.y.size(), 0.0);
    for (int i = 0; i < b.x.size(); i++)
    { 
        b_transformed.x[i] = b.x[i]/Pars::Ex;
        b_transformed.y[i] = b.y[i]/Pars::Ey;
    }

    _minBodyCoordinate[0] = *std::min_element(b_transformed.x.begin(), b_transformed.x.end()) - max_grid_size * std::sqrt(2);
    _minBodyCoordinate[1] = *std::min_element(b_transformed.y.begin(), b_transformed.y.end()) - max_grid_size * std::sqrt(2);
    _maxBodyCoordinate[0] = *std::max_element(b_transformed.x.begin(), b_transformed.x.end()) + max_grid_size * std::sqrt(2);
    _maxBodyCoordinate[1] = *std::max_element(b_transformed.y.begin(), b_transformed.y.end()) + max_grid_size * std::sqrt(2);
    
    // Determine the center and normal panel of the body
    for (int i = 0; i < _nPanel; i++)
    {
        // determine midpoint of panel
        _centerPanelX[i] = b_transformed.x[i] + ((b_transformed.x[i + 1] - b_transformed.x[i]) * 0.5e0);
        _centerPanelY[i] = b_transformed.y[i] + ((b_transformed.y[i + 1] - b_transformed.y[i]) * 0.5e0);
        // determine normal panel
        double _dx = b_transformed.x[i + 1] - b_transformed.x[i];
        double _dy = b_transformed.y[i + 1] - b_transformed.y[i];
        double _theta = atan2(_dy, _dx);
        _normalPanelX[i] = sin(_theta);  // unit normal vector in x-direction
        _normalPanelY[i] = -cos(_theta); // unit normal vector in y-direction
    }
    
    // Calculate the disctance from body to update the levels
    for (int i = 0; i < _num; i++){
        // Consider the only grid near the geometry
        double _xp = this->GridBase.x[i];
        double _yp = this->GridBase.y[i];

        if (_xp >= _minBodyCoordinate[0] && _xp <= _maxBodyCoordinate[0] &&
            _yp >= _minBodyCoordinate[1] && _yp <= _maxBodyCoordinate[1])
        {
            // initial guess value
            int _k1 = 0;
            double _rmins = 100.0e0;
            double _signNormal = 0.0e0;
            double _testPin = 0.0e0;

            for (int k = 0; k < _nPanel; k++)
            {
                // Distance between center panel and the point
                double _r2panel = std::sqrt(std::pow((_xp - _centerPanelX[k]), 2) +
                                            std::pow((_yp - _centerPanelY[k]), 2));
                if (std::abs(_r2panel) <= _rmins)
                {
                    _rmins = _r2panel; // find the center panel closest to the point
                    _k1 = k;           // Name of closest panel
                }
            }

            // normal distance from considered point to panel
            // * = 0 when grid node lie on the panel, or extended-panel
            // Negative sign means outside, positive sign mean inside
            double _rNormal = (((_xp - _centerPanelX[_k1]) * (_normalPanelX[_k1])) +
                            ((_yp - _centerPanelY[_k1]) * (_normalPanelY[_k1])));

            if (_rNormal > - max_grid_size * std::sqrt(2) && _rNormal < max_grid_size * std::sqrt(2))
                this->GridBase.level[i] = this->level_max;
        }
    }

    // Perform neighbor update
    this->level_update_neighbor();

    // Perform multiresolution adaption
    this->generate_multiblock(p);

}
*/

// [DONE] Multiresolution body adjusted distribution
void initialization::init_2d_multi_res_body_adjusted(const Body &b, Particle &par){
    
    // internal temporary variable
    double _y, _x;
    int i, j, count;
    const double baseSize = std::pow(2,Pars::max_level) * Pars::sigma;

    // ====== PARTICLE DISCRETIZATION GENERATION ====== //
    // number of particle in x and y direction
    int _nx = std::ceil(Pars::lxdom / baseSize);     // Number of particle in x direction
    int _ny = std::ceil(Pars::lydom / baseSize);     // Number of particle in y direction
    int _nparticle = _nx * _ny;                 // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);
    par.s.resize(_nparticle, baseSize);
    par.level.resize(_nparticle, 0);
    par.num = _nparticle;

    // Generate the particle position
    count = 0;
    for (i = 0; i < _nx; i++)
    {
        for (j = 0; j < _ny; j++)
        {
            _x = -Pars::xdom + (i + 0.5)*baseSize;
            _y = -baseSize*_ny*0.5 + (j + 0.5)*baseSize;
            par.x[count] = _x;
            par.y[count] = _y;
            count ++;
        }
    }

    // ====== DIVIDE AND CONQUEROR PROCEDURE ====== //
    // Perform the devide and conqueror if the point inside the finer resolution block
    int* DnC_x = new int[4]{1,-1,-1,1};   // Divide and conqueror sum for x
    int* DnC_y = new int[4]{1,1,-1,-1};   // Divide and conqueror sum for y
    
    std::vector<int> par_list;
    std::vector<int> par_list_temp;
    for (i = 0; i < par.num; i++){
        par_list.emplace_back(i);
    }

    // Constant parameter for DnC evaluation
    std::vector<double> pos;
    double R_marg;
    double R_eval;

    for (count = 0; count < Pars::max_level; count++){
        // Arithmatic Series [Sn = n/2*(2*a + b*(n-1))]
        int n_S = Pars::max_level - count;  // The position order
        double a_Sn = 15 * Pars::sigma;//0.2;//15 * Pars::sigma;     // Initial unit (Manual Adjusted)
        double b_Sn = 8 * Pars::sigma;//0.15;//8 * Pars::sigma;      // Increment (Manual adjuster)
        R_marg = n_S * (a_Sn + 0.5 * b_Sn * (n_S - 1)) * Pars::Df;

        for (auto ID:par_list){
            pos = {par.x[ID], par.y[ID]};
            R_eval = d_geom.distance_calc(pos, b);

            if (std::abs(R_eval) < R_marg)
            // DIVIDE AND CONQUEROR
            {
                // The particle is inside the finer resolution block
                _x = par.x[ID];
                _y = par.y[ID];
                
                // Point 1
                par.x[ID] = _x + 0.25 * baseSize/std::pow(2,count);
                par.y[ID] = _y + 0.25 * baseSize/std::pow(2,count);
                par.s[ID] = par.s[ID]/2;
                par.level[ID] = par.level[ID] + 1;
                par_list_temp.emplace_back(ID);

                // Point 2 - 4
                for (j = 1; j < 4; j++){
                    par_list_temp.emplace_back(par.s.size());
                    par.x.emplace_back(_x + 0.25 * baseSize/std::pow(2,count) * DnC_x[j]);
                    par.y.emplace_back(_y + 0.25 * baseSize/std::pow(2,count) * DnC_y[j]);
                    par.s.emplace_back(par.s[ID]);
                    par.level.emplace_back(par.level[ID]);
                    _nparticle ++;
                }
            }
        }
        
        // Update the particle list
        par_list = par_list_temp;
        par_list_temp.clear();
    }
    
    // Deallocate the data
    delete DnC_x, DnC_y;

    // Assign other particle properties
    par.num = _nparticle;
}

/*
Note for 2D initialization
*/

// Additional code for multibloc distribution (See at the upgraded part)
/*
void initialization::generate_multiblock(Particle &particle)
{
    // Initial the variable for calculation
    double _xp;
    double _yp;
    double _sp;
    int _lvlp;

    // The variable to store the final particle position
    std::vector<double> _xB;
    std::vector<double> _yB;
    std::vector<double> _sB;
    std::vector<int> _lvlB;
    std::vector<int> _lblB;
    
    // bool _isFinish = false; // Initiate

    for (size_t i = 0; i < this->GridBase.x.size(); i++)
    {
        // Assign the variable
        _xp = this->GridBase.x[i];
        _yp = this->GridBase.y[i];
        _sp = this->GridBase.s[i];
        _lvlp = this->GridBase.level[i];

        // Determine the number of iteration needed
        int iter_num = _lvlp - 1;

        // The variable to store the temporarry particle position for sub level
        std::vector<double> _xBcalc(1,_xp);
        std::vector<double> _yBcalc(1,_yp);

        // The variable to store the temporarry particle position for calculation
        std::vector<double> _xBtemp;
        std::vector<double> _yBtemp;

        // Iteration of devide and conqueror for each position
        for (int iter = 0; iter < iter_num; iter ++)
        {
            _xBtemp.clear();
            _yBtemp.clear();
            
            for (int k = 0; k < _xBcalc.size(); k++)
            {
                _xBtemp.push_back(_xBcalc[k] + 0.25 * _sp);
                _xBtemp.push_back(_xBcalc[k] + 0.25 * _sp);
                _xBtemp.push_back(_xBcalc[k] - 0.25 * _sp);
                _xBtemp.push_back(_xBcalc[k] - 0.25 * _sp);

                _yBtemp.push_back(_yBcalc[k] + 0.25 * _sp);
                _yBtemp.push_back(_yBcalc[k] - 0.25 * _sp);
                _yBtemp.push_back(_yBcalc[k] + 0.25 * _sp);
                _yBtemp.push_back(_yBcalc[k] - 0.25 * _sp);
            }
            _xBcalc = _xBtemp;
            _yBcalc = _yBtemp;
            _sp = 0.5 * _sp;
        }
        
        // Store the final calculated particle position to the particle storage
        for (int k = 0; k < _xBcalc.size(); k++){
            _xB.push_back(_xBcalc[k]);
            _yB.push_back(_yBcalc[k]);
            _sB.push_back(_sp);
            _lvlB.push_back(_lvlp);
            _lblB.push_back(i);
        }
    }

    // Still not appropriated
    particle.x = _xB;
    particle.y = _yB;
    particle.s = _sB;
    particle.level = _lvlB;
    particle.label = _lblB;
    particle.num = _xB.size();
    particle.u.resize(particle.num, 0.0e0);
    particle.v.resize(particle.num, 0.0e0);
    particle.gz.resize(particle.num, 0.0e0);
    particle.vorticity.resize(particle.num, 0.0e0);
    particle.neighbor.resize(particle.num, std::vector<int>());
    particle.isActive.resize(particle.num, false);
    particle.chi.resize(particle.num, 0.0);
}
*/
