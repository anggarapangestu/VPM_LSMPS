#include "adaptation.hpp"
#include "../grid_block/gridNodeNgh.hpp"
#include "../LSMPS/LSMPSa.hpp"
#include "../save_data/save_data.hpp"

#define ERROR_ADJUSTMENT false       // A selection to adjust the error profile
#define ERROR_ADJUSTMENT_TYPE 2     // Adjusting the error profile with
                                    //  1:= Direct flat profile
                                    //  2:= A smoothing profile

/**
 *  @brief  Particle adaptation manager distribution manager. Generate the particle data
 *  using the selected type of distribution by user in global.cpp file.
 *         
 *  @param  _particle  Particle data contains the physical properties.
 *  @param  _bodyList  The list of body data container used for particle generation.
 *  @param  _gridNode  Grid node data container used for particle generation.
*/
void adaptation::get_adaptation_LSMPS(Particle &_parEval, GridNode &_baseGrid){
    // Perform particle adaptation
    // ***************************

    // PATCH UPDATE:
    // PROCEDURE 1:
    // ***********
    // Calculate the properties selection for the error predictor
    this->predProp.clear();
    // this->error_selector(_parEval);

    // Create the error properties selector (Use a global container declaration)
    std::vector<double>* lap;          // The pointer for vorticity value
    lap = new std::vector<double>;     // Create a dynamic memory for saving

    // Assign the value into the container
    lap->clear(); lap->resize(_parEval.num);

    if (Pars::opt_adapt_err_pred == 1)
    {
        #pragma omp parallel for
        for (int i = 0; i < _parEval.num; i++){
            (*lap)[i] = std::abs(_parEval.F_a[i]);
        }
    }else if(Pars::opt_adapt_err_pred == 2){
        #pragma omp parallel for
        for (int i = 0; i < _parEval.num; i++){
            (*lap)[i] = std::abs(_parEval.Fx2_a[i]) + std::abs(_parEval.Fy2_a[i]);
        }
    }else if(Pars::opt_adapt_err_pred == 3){
        #pragma omp parallel for
        for (int i = 0; i < _parEval.num; i++){
            (*lap)[i] = std::abs(_parEval.Fx2_a[i] + _parEval.Fy2_a[i]);
        }
    }


    // Done the calculation, update the container value
    this->predProp.push_back(lap);

    // PROCEDURE 2:
    // ***********
    // Perform the Particle Adaptation using gridNode
    GridNodeAdapt _gridNodeTools;
    bool isAdapted;
    _gridNodeTools.get_adaptation_LSMPS(_baseGrid, _parEval, this->predProp);

    // Free the predictor properties
    for (int i = 0; i < this->predProp.size(); i++){
        delete this->predProp[i];
    }

    return;
}

/**
 *  @brief  Particle adaptation manager distribution manager. Generate the particle data
 *  using the selected type of distribution by user in global.cpp file.
 *         
 *  @param  _particle  Particle data contains the physical properties.
 *  @param  _bodyList  The list of body data container used for particle generation.
 *  @param  _gridNode  Grid node data container used for particle generation.
*/
bool adaptation::get_adaptation(Particle &_parEval, Particle *&_parBase, GridNode &_baseGrid){
    // Perform particle adaptation
    // ***************************

    // // PATCH UPDATE:
    // // PROCEDURE 1:
    // // ***********
    // // First assign the lagrangian particle into the current base grid
    // std::unordered_map<int, std::vector<int>> parNodeMap;       // The particle node ID map (the contract between lagrange particle to eulerian grid)
    // GridNodeNgh grdNghTool;
    // grdNghTool.assign_par2node(_baseGrid, parNodeMap, _parEval);

    
    // PROCEDURE 2:
    // ***********
    // Calculate the properties selection for the error predictor
    this->predProp.clear();
    this->error_selector(_parEval);

    // std::ofstream writter;
    // writter.open("output/errSel.csv");
    // #if (DIM == 2)
    // writter << "x,y,s";
    // #elif (DIM == 3)
    // writter << "x,y,z,s";
    // #endif
    // for (int i = 0 ; i < this->predProp.size(); i++){
    //     writter << ",p" << i+1;
    // }
    // writter << "\n";
    // for (int i = 0; i < _parEval.num; i++){
    //     writter << _parEval.x[i] << ","
    //             << _parEval.y[i] << ",";
    //     #if (DIM == 3)
    //     writter << _parEval.z[i] << ",";
    //     #endif
    //     writter << _parEval.s[i];
    //     for (int k = 0 ; k < this->predProp.size(); k++){
    //         writter << "," << (*this->predProp[k])[i];
    //     }
    //     writter << "\n";
    // }
    // writter.close();
    // exit(0);


    // PROCEDURE 3:
    // ***********
    // Perform the Particle Adaptation using gridNode
    GridNodeAdapt _gridNodeTools;
    bool isAdapted;
    // isAdapted = _gridNodeTools.get_adaptation(_baseGrid, *_parBase, _parEval, 0, this->predProp, parNodeMap);
    isAdapted = _gridNodeTools.get_adaptation(_baseGrid, *_parBase, _parEval, this->predProp);

    // _baseGrid.saveLeafGrid(_baseGrid, "After_Adapted");

    // Free the predictor properties
    for (int i = 0; i < this->predProp.size(); i++){
        delete this->predProp[i];
    }
    
    // Check if there is no adaptation
    if (isAdapted == false){
        return isAdapted;
    }
    
    // Retrieve the data to base particle
    delete _parBase;        // Avoid memory leaking (The old data template particle is not necessary)
    _gridNodeTools.take_particle_pointer(_parBase);

    // Update the particle data
    // ************************
    int &num = _parBase->num;
    // Resize the other particle data
    // [NOTE] Modify if there is any properties update
    _parBase->u.resize(num);
    _parBase->v.resize(num);
    if (DIM == 3){
        _parBase->w.resize(num);
        _parBase->vortx.resize(num);
        _parBase->vorty.resize(num);
        _parBase->vortz.resize(num);
    }
    _parBase->gz.resize(num);
    _parBase->vorticity.resize(num);

    return isAdapted;
}

/**
 *  @brief  Particle adaptation error predictor type selector. This is a manager function.
 *  Use to select the type for adaptation error predictor.
 *         
 *  @param  _parEval  Particle data container.
 *  @param  _type  The type of error predictor selected [1:= Feature Predictor; 2:= Gradient Predictor]
*/
void adaptation::error_selector(const Particle &_parEval){
    // Evaluate the properties for error predictor evaluation
    switch (Pars::opt_adapt_err_pred){
    case 1:	// Feature base
        printf("%sType 1: Feature Based on Vorticity %s\n", FONT_CYAN, FONT_RESET);
        this->feature_predictor(_parEval);
        break;
    case 2:	// Gradient base
        printf("%sType 2: Gradient Based on Vorticity %s\n", FONT_CYAN, FONT_RESET);
        this->gradient_predictor(_parEval);
        break;
    case 3:	// Laplacian base
        printf("%sType 3: Laplacian Based for LSMPS %s\n", FONT_CYAN, FONT_RESET);
        this->laplace_predictor(_parEval);
        break;
    case 4: 
        // Reserved ... 
    default:
        break;
    }

    // Adjust the head profile of the error
    #if (ERROR_ADJUSTMENT)
        this->error_adjustment(0);
    #endif

    return;
}

// Adjust the error selector
void adaptation::error_adjustment(int winPos){
    // Adjust all error profile at all winPos
    /**
     * The configuration of tolerance windows (By only changing the window 0)
     *   =============================== Max            =============================== Max    
     *      Window 0                                       Window 0                            
     *   --------------------------- n^1 Max                             Expanded              
     *      Window 1                           __|.     --------------------------- n^1 Max*
     *   --------------------------- n^2 Max  |__  >       Window 1                        
     *      Window 2                             |'     --------------------------- n^2 Max*
     *   --------------------------- n^2 Max               Window 2                        
     *      Window 3                                    --------------------------- n^2 Max*
     *   --------------------------- n^3 Max                ... etc.
     *     ... etc.
    */

    
    // Intermediate variable
    int numPred = this->predProp.size();
    std::vector<double> maxVal(numPred);

    // Evaluate each maximum value
    for (int i = 0; i < numPred; i++){
        maxVal[i] = *std::max_element(this->predProp[i]->begin(),this->predProp[i]->end());
    }

    // A detailed check
    std::vector<bool> isWork(numPred, true);
    for (int i = 0; i < numPred; i++){
        if (maxVal[i] < __DBL_EPSILON__){
            isWork[i] = false;
        }
    }


    // Evaluation on the error profile
    #if (ERROR_ADJUSTMENT_TYPE == 1)
        // For current convenient we can just drop it to the given tolerance
        for (int i = 0; i < numPred; i++){
            // Alias to the current property
            std::vector<double>& _prop = *this->predProp.at(i);
            double _maxVal = maxVal[i];

            // Calcualte the current properties
            for (auto& val : _prop){
                double valRel = val/_maxVal;

                // Check whether the current 
                if (valRel > Pars::adapt_head_tol){
                    val = _maxVal * Pars::adapt_head_tol;
                }
            }
        }
    #elif (ERROR_ADJUSTMENT_TYPE == 2)
        // The smoothing stretch of error profile
        
        // The common variable in the following calculation
        double n_h = Pars::adapt_head_tol;
        double n_t = Pars::adapt_tol;
        int k = winPos;
        double tolBorder = n_h * std::pow(n_t, k);

        // Evalaute all error properties
        for (int i = 0; i < numPred; i++){
            // Alias to the current property
            std::vector<double>& _prop = *this->predProp.at(i);
            double _maxVal = maxVal[i];

            // Calcualte the current properties
            for (auto& val : _prop){
                double valRel = val/_maxVal;

                // Check whether the current 
                if (valRel > tolBorder){
                    double _grad = std::log10(std::pow(n_t,k+1)) / (std::log10(n_h * std::pow(n_t,k)));
                    double _X = std::log10(valRel);
                    double _B = std::log10(n_h/n_t);
                    double _exp = _grad*_X + _B;
                    val = _maxVal * std::pow(10.0, _exp);
                }
            }
        }
    #endif

    return;
}


/**
 *  @brief  The particle evaluation for [feature type] error predictor adaptation.
 *         
 *  @param  _parEval  Particle data container.
*/
void adaptation::feature_predictor(const Particle &_parEval){
    // Evaluate the properties for error predictor evaluation

    // Create the error properties selector (Use a global container declaration)
    std::vector<double>* vort;          // The pointer for vorticity value
    vort = new std::vector<double>;     // Create a dynamic memory for saving

    // Assign the value into the container
    vort->clear(); vort->resize(_parEval.num);

    #pragma omp parallel for
    for (int i = 0; i < _parEval.num; i++){
        (*vort)[i] = std::abs(_parEval.vorticity[i]);
    }

    // Done the calculation, update the container value
    this->predProp.push_back(vort);

    return;
}

/**
 *  @brief  The particle evaluation for [feature type] error predictor adaptation.
 *         
 *  @param  _parEval  Particle data container.
*/
void adaptation::laplace_predictor(const Particle &_parEval){
    // Evaluate the properties for error predictor evaluation

    // A baypass Method
    #if DIM == 2  // Calculate the gradient for 2 Dimension
        if(!_parEval.absLapVortZ.empty()){
            // Calculate the bypass method then completed

            // Create the error properties selector (Use a global container declaration)
            std::vector<double>* gradVort;          // The pointer for vorticity gradient value
            gradVort = new std::vector<double>;     // Create a dynamic memory for saving

            // Assign the value into the container
            gradVort->clear(); gradVort->resize(_parEval.num, 0.0);
            #pragma omp parallel for
            for (int i = 0; i < _parEval.num; i++){
                (*gradVort)[i] = _parEval.absLapVortZ[i];
            }

            // Done the calculation, update the container value
            this->predProp.push_back(gradVort);

            std::cout << "Error Properties Calculation is Done by BYPASS!!\n";
            // double _A;
            // std::cin >> _A;
            return;
        }
    #elif DIM == 3    // Calculate the gradient for 3 Dimension
        if( !_parEval.absLapVortX.empty() && 
            !_parEval.absLapVortY.empty() && 
            !_parEval.absLapVortZ.empty() ){
            
            // Calculate the bypass method then completed

            // Create the error properties selector (Use a global container declaration)
            std::vector<double>* gradVortX;          // The pointer for vorticity gradient value
            std::vector<double>* gradVortY;          // The pointer for vorticity gradient value
            std::vector<double>* gradVortZ;          // The pointer for vorticity gradient value
            gradVortX = new std::vector<double>;     // Create a dynamic memory for saving
            gradVortY = new std::vector<double>;     // Create a dynamic memory for saving
            gradVortZ = new std::vector<double>;     // Create a dynamic memory for saving

            // Assign the value into the container
            gradVortX->clear(); gradVortX->resize(_parEval.num, 0.0);
            gradVortY->clear(); gradVortY->resize(_parEval.num, 0.0);
            gradVortZ->clear(); gradVortZ->resize(_parEval.num, 0.0);
            #pragma omp parallel for
            for (int i = 0; i < _parEval.num; i++){
                (*gradVortX)[i] = _parEval.absLapVortX[i];
                (*gradVortY)[i] = _parEval.absLapVortY[i];
                (*gradVortZ)[i] = _parEval.absLapVortZ[i];
            }

            // // Done the calculation, update the container value
            // this->predProp.push_back(gradVortX);
            // this->predProp.push_back(gradVortY);
            // this->predProp.push_back(gradVortZ);

            // Done the calculation, update the container value
            MESSAGE_LOG << "Error Predictor Evaluation!\n";
            if (*std::max_element(gradVortX->begin(),gradVortX->end()) > __DBL_EPSILON__){
                this->predProp.push_back(gradVortX);
            }else{
                std::cout << "Laplacian of VortX is not included!\n";
            }
            
            if (*std::max_element(gradVortY->begin(),gradVortY->end()) > __DBL_EPSILON__){
                this->predProp.push_back(gradVortY);
            }else{
                std::cout << "Laplacian of VortY is not included!\n";
            }

            if (*std::max_element(gradVortZ->begin(),gradVortZ->end()) > __DBL_EPSILON__){
                this->predProp.push_back(gradVortZ);
            }else{
                std::cout << "Laplacian of VortZ is not included!\n";
            }

            // std::cout << "Error Properties Calculation is Done by BYPASS!!\n";
            // double _A;
            // std::cin >> _A;
            return;
        }

    #endif

    // Recalculation method (Too much calculation)
    
    // PROCEDURE 1:
    // ************
    // Collecting the only active particle data
    
    // Take the active particle only
    // Internal variables
    const Particle& p = _parEval;   // The data of particle inside 'active particle + buffer zone'
    Particle _particle;             // The data of particle inside 'active particle + buffer zone'
    std::vector<int> _index;        // The index list of particle inside 'active particle + buffer zone'
    std::vector<bool> _flag(p.num, false);      // The flag list of evaluated particle, At initial still no evaluated particle list

    // Assign the index of active particles
    for (int ID = 0; ID < p.num; ID++)
    {
        // Check if the particle ID is active
        if (p.isActive[ID] == true)
        {
            // [1] The current ID particle
            if(!Pars::flag_ngh_include_self && _flag[ID] == false)
            {
                // Assign ID into index list if not in the list (if flag == false)
                _index.push_back(ID);    // Assign ID into the index list
                _flag[ID] = true;        // Update the evaluation flag of particle ID
            }

            // [2] The neighbor particles of current ID particle
            for (auto _ngh_ID : p.neighbor[ID])
            {
                if(_flag[_ngh_ID] == false)
                {
                    // Assign into index list if not in the list (if flag == false)
                    _index.push_back(_ngh_ID);  // Assign _ngh_ID particle into the index list
                    _flag[_ngh_ID] = true;      // Update the evaluation flag of particle _ngh_ID
                }
            }
        }
    }

    // Resize the _particle variable
    _particle.num = _index.size();
    _particle.x.clear(); _particle.x.resize(_particle.num);
    _particle.y.clear(); _particle.y.resize(_particle.num);
    #if DIM == 2
    _particle.vorticity.clear(); _particle.vorticity.resize(_particle.num);
    #elif DIM == 3
    _particle.z.clear(); _particle.z.resize(_particle.num);
    _particle.vortx.clear(); _particle.vortx.resize(_particle.num);
    _particle.vorty.clear(); _particle.vorty.resize(_particle.num);
    _particle.vortz.clear(); _particle.vortz.resize(_particle.num);
    #endif
    _particle.s.clear(); _particle.s.resize(_particle.num);
    _particle.neighbor.clear(); _particle.neighbor.resize(_particle.num);

    // Store the particle data of each particle inside _index list
    #pragma omp parallel for
    for (int i = 0; i < _particle.num; i++)
    {
        const int &_ID = _index[i];
        _particle.s[i] = (p.s[_ID]);
        _particle.x[i] = (p.x[_ID]);
        _particle.y[i] = (p.y[_ID]);
        #if DIM == 2
        _particle.vorticity[i] = (p.vorticity[_ID]);
        #elif DIM == 3
        _particle.z[i] = (p.z[_ID]);
        _particle.vortx[i] = (p.vortx[_ID]);
        _particle.vorty[i] = (p.vorty[_ID]);
        _particle.vortz[i] = (p.vortz[_ID]);
        #endif
        _particle.neighbor[i] = (p.neighbor[_ID]);
    }

    
    // PROCEDURE 2:
    // ***********
    // Calculate the gradient
    #if DIM == 2  // Calculate the gradient for 2 Dimension
        // Create the error properties selector (Use a global container declaration)
        std::vector<double>* gradVort;          // The pointer for vorticity gradient value
        gradVort = new std::vector<double>;     // Create a dynamic memory for saving

        // Internal method
        LSMPSa lsmpsa;		          // To calculate laplacian

        // Calculate the second order differential of vorticity
        Particle &_p = _particle;       // Aliasing
        lsmpsa.set_LSMPS(_p.x, _p.y, _p.s, _p.vorticity, p.x, p.y, p.s, p.vorticity, _p.neighbor);
        std::vector<double> _d2fd2x = lsmpsa.get_d2d2x();
        std::vector<double> _d2fd2y = lsmpsa.get_d2d2y();

        // Assign the value into the container
        gradVort->clear(); gradVort->resize(_parEval.num, 0.0);
        #pragma omp parallel for
        for (int i = 0; i < /*_parEval.num*/_index.size(); i++){
            const int &_ID = _index[i];
            (*gradVort)[_ID] = std::abs(_d2fd2x[i] + _d2fd2y[i]);
        }

        // Done the calculation, update the container value
        this->predProp.push_back(gradVort);
        
    #elif DIM == 3    // Calculate the gradient for 3 Dimension
        // Create the error properties selector (Use a global container declaration)
        std::vector<double>* gradVortX;          // The pointer for vorticity gradient value
        std::vector<double>* gradVortY;          // The pointer for vorticity gradient value
        std::vector<double>* gradVortZ;          // The pointer for vorticity gradient value
        gradVortX = new std::vector<double>;     // Create a dynamic memory for saving
        gradVortY = new std::vector<double>;     // Create a dynamic memory for saving
        gradVortZ = new std::vector<double>;     // Create a dynamic memory for saving

        // Internal method
        LSMPSa lsmpsa;		          // To calculate laplacian

        // Calculate the second order differential of vorticity
        Particle &_p = _particle;       // Aliasing
        // Calculate the x directing vorticity
        lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vortx, 
                            p.x, p.y, p.z, p.s, p.vortx, _p.neighbor);
        std::vector<double> _d2vortxd2x = lsmpsa.get_d2d2x();
        std::vector<double> _d2vortxd2y = lsmpsa.get_d2d2y();
        std::vector<double> _d2vortxd2z = lsmpsa.get_d2d2z();

        // Calculate the y directing vorticity
        lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vorty, 
                            p.x, p.y, p.z, p.s, p.vorty, _p.neighbor);
        std::vector<double> _d2vortyd2x = lsmpsa.get_d2d2x();
        std::vector<double> _d2vortyd2y = lsmpsa.get_d2d2y();
        std::vector<double> _d2vortyd2z = lsmpsa.get_d2d2z();

        // Calculate the z directing vorticity
        lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vortz, 
                            p.x, p.y, p.z, p.s, p.vortz, _p.neighbor);
        std::vector<double> _d2vortzd2x = lsmpsa.get_d2d2x();
        std::vector<double> _d2vortzd2y = lsmpsa.get_d2d2y();
        std::vector<double> _d2vortzd2z = lsmpsa.get_d2d2z();

        // Assign the value into the container
        gradVortX->clear(); gradVortX->resize(_parEval.num, 0.0);
        gradVortY->clear(); gradVortY->resize(_parEval.num, 0.0);
        gradVortZ->clear(); gradVortZ->resize(_parEval.num, 0.0);
        #pragma omp parallel for
        for (int i = 0; i < /*_parEval.num*/_index.size(); i++){
            const int &_ID = _index[i];
            (*gradVortX)[_ID] = std::abs(_d2vortxd2x[i] + _d2vortxd2y[i] + _d2vortxd2z[i]);
            (*gradVortY)[_ID] = std::abs(_d2vortyd2x[i] + _d2vortyd2y[i] + _d2vortyd2z[i]);
            (*gradVortZ)[_ID] = std::abs(_d2vortzd2x[i] + _d2vortzd2y[i] + _d2vortzd2z[i]);
        }

        // Done the calculation, update the container value
        MESSAGE_LOG << "Error Predictor Evaluation!\n";
        if (*std::max_element(gradVortX->begin(),gradVortX->end()) > __DBL_EPSILON__){
            this->predProp.push_back(gradVortX);
        }else{
            std::cout << "Laplacian of VortX is not included!\n";
        }
        
        if (*std::max_element(gradVortY->begin(),gradVortY->end()) > __DBL_EPSILON__){
            this->predProp.push_back(gradVortY);
        }else{
            std::cout << "Laplacian of VortY is not included!\n";
        }

        if (*std::max_element(gradVortZ->begin(),gradVortZ->end()) > __DBL_EPSILON__){
            this->predProp.push_back(gradVortZ);
        }else{
            std::cout << "Laplacian of VortZ is not included!\n";
        }
    #endif

    // // Done the calculation
    std::cout << "Error Properties Calculation is Done by RECALCULATION!!\n";
    // double _A;
    // std::cin >> _A;
    return;
}

/**
 *  @brief  The particle evaluation for [gradient type] error predictor adaptation.
 *  The gradient using the 2nd second order gradient
 *         
 *  @param  _parEval  Particle data container of physical properties.
*/
void adaptation::gradient_predictor(const Particle &_parEval){
    // Evaluate the properties for error predictor evaluation

    // A baypass Method
    #if DIM == 2  // Calculate the gradient for 2 Dimension
        if(!_parEval.absLapVortZ.empty()){
            // Calculate the bypass method then completed

            // Create the error properties selector (Use a global container declaration)
            std::vector<double>* gradVort;          // The pointer for vorticity gradient value
            gradVort = new std::vector<double>;     // Create a dynamic memory for saving

            // Assign the value into the container
            gradVort->clear(); gradVort->resize(_parEval.num, 0.0);
            #pragma omp parallel for
            for (int i = 0; i < _parEval.num; i++){
                (*gradVort)[i] = _parEval.absLapVortZ[i];
            }

            // Done the calculation, update the container value
            this->predProp.push_back(gradVort);

            // std::cout << "Error Properties Calculation is Done by BYPASS!!\n";
            // double _A;
            // std::cin >> _A;
            return;
        }
    #elif DIM == 3    // Calculate the gradient for 3 Dimension
        if( !_parEval.absLapVortX.empty() && 
            !_parEval.absLapVortY.empty() && 
            !_parEval.absLapVortZ.empty() ){
            
            // Calculate the bypass method then completed

            // Create the error properties selector (Use a global container declaration)
            std::vector<double>* gradVortX;          // The pointer for vorticity gradient value
            std::vector<double>* gradVortY;          // The pointer for vorticity gradient value
            std::vector<double>* gradVortZ;          // The pointer for vorticity gradient value
            gradVortX = new std::vector<double>;     // Create a dynamic memory for saving
            gradVortY = new std::vector<double>;     // Create a dynamic memory for saving
            gradVortZ = new std::vector<double>;     // Create a dynamic memory for saving

            // Assign the value into the container
            gradVortX->clear(); gradVortX->resize(_parEval.num, 0.0);
            gradVortY->clear(); gradVortY->resize(_parEval.num, 0.0);
            gradVortZ->clear(); gradVortZ->resize(_parEval.num, 0.0);
            #pragma omp parallel for
            for (int i = 0; i < _parEval.num; i++){
                (*gradVortX)[i] = _parEval.absLapVortX[i];
                (*gradVortY)[i] = _parEval.absLapVortY[i];
                (*gradVortZ)[i] = _parEval.absLapVortZ[i];
            }

            // Done the calculation, update the container value
            this->predProp.push_back(gradVortX);
            this->predProp.push_back(gradVortY);
            this->predProp.push_back(gradVortZ);

            // std::cout << "Error Properties Calculation is Done by BYPASS!!\n";
            // double _A;
            // std::cin >> _A;
            return;
        }

    #endif

    // Recalculation method (Too much calculation)
    
    // PROCEDURE 1:
    // ************
    // Collecting the only active particle data
    
    // Take the active particle only
    // Internal variables
    const Particle& p = _parEval;   // The data of particle inside 'active particle + buffer zone'
    Particle _particle;             // The data of particle inside 'active particle + buffer zone'
    std::vector<int> _index;        // The index list of particle inside 'active particle + buffer zone'
    std::vector<bool> _flag(p.num, false);      // The flag list of evaluated particle, At initial still no evaluated particle list

    // Assign the index of active particles
    for (int ID = 0; ID < p.num; ID++)
    {
        // Check if the particle ID is active
        if (p.isActive[ID] == true)
        {
            // [1] The current ID particle
            if(!Pars::flag_ngh_include_self && _flag[ID] == false)
            {
                // Assign ID into index list if not in the list (if flag == false)
                _index.push_back(ID);    // Assign ID into the index list
                _flag[ID] = true;        // Update the evaluation flag of particle ID
            }

            // [2] The neighbor particles of current ID particle
            for (auto _ngh_ID : p.neighbor[ID])
            {
                if(_flag[_ngh_ID] == false)
                {
                    // Assign into index list if not in the list (if flag == false)
                    _index.push_back(_ngh_ID);  // Assign _ngh_ID particle into the index list
                    _flag[_ngh_ID] = true;      // Update the evaluation flag of particle _ngh_ID
                }
            }
        }
    }

    // Resize the _particle variable
    _particle.num = _index.size();
    _particle.x.clear(); _particle.x.resize(_particle.num);
    _particle.y.clear(); _particle.y.resize(_particle.num);
    #if DIM == 2
    _particle.vorticity.clear(); _particle.vorticity.resize(_particle.num);
    #elif DIM == 3
    _particle.z.clear(); _particle.z.resize(_particle.num);
    _particle.vortx.clear(); _particle.vortx.resize(_particle.num);
    _particle.vorty.clear(); _particle.vorty.resize(_particle.num);
    _particle.vortz.clear(); _particle.vortz.resize(_particle.num);
    #endif
    _particle.s.clear(); _particle.s.resize(_particle.num);
    _particle.neighbor.clear(); _particle.neighbor.resize(_particle.num);

    // Store the particle data of each particle inside _index list
    #pragma omp parallel for
    for (int i = 0; i < _particle.num; i++)
    {
        const int &_ID = _index[i];
        _particle.s[i] = (p.s[_ID]);
        _particle.x[i] = (p.x[_ID]);
        _particle.y[i] = (p.y[_ID]);
        #if DIM == 2
        _particle.vorticity[i] = (p.vorticity[_ID]);
        #elif DIM == 3
        _particle.z[i] = (p.z[_ID]);
        _particle.vortx[i] = (p.vortx[_ID]);
        _particle.vorty[i] = (p.vorty[_ID]);
        _particle.vortz[i] = (p.vortz[_ID]);
        #endif
        _particle.neighbor[i] = (p.neighbor[_ID]);
    }

    
    // PROCEDURE 2:
    // ***********
    // Calculate the gradient
    #if DIM == 2  // Calculate the gradient for 2 Dimension
        // Create the error properties selector (Use a global container declaration)
        std::vector<double>* gradVort;          // The pointer for vorticity gradient value
        gradVort = new std::vector<double>;     // Create a dynamic memory for saving

        // Internal method
        LSMPSa lsmpsa;		          // To calculate laplacian

        // Calculate the second order differential of vorticity
        Particle &_p = _particle;       // Aliasing
        lsmpsa.set_LSMPS(_p.x, _p.y, _p.s, _p.vorticity, p.x, p.y, p.s, p.vorticity, _p.neighbor);
        std::vector<double> _d2fd2x = lsmpsa.get_d2d2x();
        std::vector<double> _d2fd2y = lsmpsa.get_d2d2y();

        // Assign the value into the container
        gradVort->clear(); gradVort->resize(_parEval.num, 0.0);
        #pragma omp parallel for
        for (int i = 0; i < /*_parEval.num*/_index.size(); i++){
            const int &_ID = _index[i];
            (*gradVort)[_ID] = std::abs(_d2fd2x[i]) + std::abs(_d2fd2y[i]);
        }

        // Done the calculation, update the container value
        this->predProp.push_back(gradVort);
        
    #elif DIM == 3    // Calculate the gradient for 3 Dimension
        // Create the error properties selector (Use a global container declaration)
        std::vector<double>* gradVortX;          // The pointer for vorticity gradient value
        std::vector<double>* gradVortY;          // The pointer for vorticity gradient value
        std::vector<double>* gradVortZ;          // The pointer for vorticity gradient value
        gradVortX = new std::vector<double>;     // Create a dynamic memory for saving
        gradVortY = new std::vector<double>;     // Create a dynamic memory for saving
        gradVortZ = new std::vector<double>;     // Create a dynamic memory for saving

        // Internal method
        LSMPSa lsmpsa;		          // To calculate laplacian

        // Calculate the second order differential of vorticity
        Particle &_p = _particle;       // Aliasing
        // Calculate the x directing vorticity
        lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vortx, 
                            p.x, p.y, p.z, p.s, p.vortx, _p.neighbor);
        std::vector<double> _d2vortxd2x = lsmpsa.get_d2d2x();
        std::vector<double> _d2vortxd2y = lsmpsa.get_d2d2y();
        std::vector<double> _d2vortxd2z = lsmpsa.get_d2d2z();

        // Calculate the y directing vorticity
        lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vorty, 
                            p.x, p.y, p.z, p.s, p.vorty, _p.neighbor);
        std::vector<double> _d2vortyd2x = lsmpsa.get_d2d2x();
        std::vector<double> _d2vortyd2y = lsmpsa.get_d2d2y();
        std::vector<double> _d2vortyd2z = lsmpsa.get_d2d2z();

        // Calculate the z directing vorticity
        lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vortz, 
                            p.x, p.y, p.z, p.s, p.vortz, _p.neighbor);
        std::vector<double> _d2vortzd2x = lsmpsa.get_d2d2x();
        std::vector<double> _d2vortzd2y = lsmpsa.get_d2d2y();
        std::vector<double> _d2vortzd2z = lsmpsa.get_d2d2z();

        // Assign the value into the container
        gradVortX->clear(); gradVortX->resize(_parEval.num, 0.0);
        gradVortY->clear(); gradVortY->resize(_parEval.num, 0.0);
        gradVortZ->clear(); gradVortZ->resize(_parEval.num, 0.0);
        #pragma omp parallel for
        for (int i = 0; i < /*_parEval.num*/_index.size(); i++){
            const int &_ID = _index[i];
            (*gradVortX)[_ID] = std::abs(_d2vortxd2x[i]) + std::abs(_d2vortxd2y[i]) + std::abs(_d2vortxd2z[i]);
            (*gradVortY)[_ID] = std::abs(_d2vortyd2x[i]) + std::abs(_d2vortyd2y[i]) + std::abs(_d2vortyd2z[i]);
            (*gradVortZ)[_ID] = std::abs(_d2vortzd2x[i]) + std::abs(_d2vortzd2y[i]) + std::abs(_d2vortzd2z[i]);
        }

        // Done the calculation, update the container value
        this->predProp.push_back(gradVortX);
        this->predProp.push_back(gradVortY);
        this->predProp.push_back(gradVortZ);
    #endif

    // // Done the calculation
    // std::cout << "Error Properties Calculation is Done by RECALCULATION!!\n";
    // double _A;
    // std::cin >> _A;
    return;
}

