#include "gridNodeAdapt.hpp"

// A flag to also adapt the resolution inside the body
#define ADAPT_INSIDE_BODY 1

// ===================================================== //
// +------------ Adaptation Initialization ------------+ //
// ===================================================== //
// #pragma region INITIALIZATION

/**
 *  @brief  Set up the target resolution level of each LEAF node.
 *         
 *  @param  baseGrid [UPDATE] The node container which node TRL is to be updated.
 *  @param  tempPar  The template particle to retrieve particle coordinate data.
 *  @param  errPredProp  The container of the selected properties of each particle, sorted by particle ID.
*/
void GridNodeAdapt::targetLevel(GridNode &baseGrid, const Particle &tempPar, const std::vector<std::vector<double>*>& errPredProp){
    // In this method, the value of "selected properties" of each particle are being compared to the MAX_VALUE
    // The target resolution level of the given particle takes the matched windows, see the following illustration

    /* ILLUSTRATION OF Target Resolution Level (TRL)
        For the sake of this illustration, take an example value for each parameter as follow.
        let. MAX_LEVEL = 3, MAX_VALUE = 5.0, TOLERANCE = 0.1, and VALUE = 0.1
        
        [LARGEST VALUE]
          VAL 0 -----------------------------> 5.000 :(MAX_VALUE * (TOLERANCE)^0)
            ^
            |    WINDOW of TRL -> [3] (LEAF) : MAX_LEVEL
            |
          VAL 1 -----------------------------> 0.500 :(MAX_VALUE * (TOLERANCE)^1)
            ^
            |    WINDOW of TRL -> [2] : MAX_LEVEL - 1
            |        [x] The particle value [0.1] lies in this window
            |     thus the TRL of current particle is [2]
            |
          VAL 2 -----------------------------> 0.050 :(MAX_VALUE * (TOLERANCE)^2)
            ^    
            |    WINDOW of TRL -> [1] : MAX_LEVEL - 2
            |
          VAL 3 -----------------------------> 0.005 :(MAX_VALUE * (TOLERANCE)^3)
            ^    
            |    WINDOW of TRL -> [0] (ROOT) : MAX_LEVEL - 3
            |
          ABS 0 -----------------------------> 0
        [LOWEST VALUE]
    */

    // // Internal variable
    // int tarLv;      // The target resolution level of evaluated particle
    // const double winTol = std::log(Pars::adapt_tol);      // The tolerance value for windows

    // Set all nodes tarResLv in the GridNode to ROOT prior to the evaluation
    for (auto &[_ID,_node] : baseGrid.nodeMap){
        // Reset the node target level
        _node->tarResLv = ROOT_LEVEL;
    }

    // // // [DEBUG LINE]
    // // std::vector<int> TARGET_LVL_LIST(this->maxLevel+1,0);
    // const double intial_tol = 5.0e-1;       // Additional patch added at 28 May 2024

    // Section 1
    // *********
    // Evaluate the maximum value of each 
    this->maxValue.resize(errPredProp.size());
    // double max_value [errPredProp.size()];
    for (int d = 0; d < errPredProp.size(); d++){
        this->maxValue[d] = *(std::max_element(errPredProp[d]->begin(), errPredProp[d]->end()));
    }
    
    // for (int i = 0; i < this->maxValue.size(); i++){
    //     std::cout << "The maximum value is " << this->maxValue[i] << "\n";
    // }
    
    // Section 2
    // *********
    // Evaluate through all particle
    for (int _parID = 0; _parID < tempPar.num; _parID++){
        
        // // Iterative Evaluation
        // // ********************
        // // Check the TRL of current particle
        // double const &currValue = selProp[_parID];
        // // double tolerance = this->maxValue * intial_tol/*Pars::adapt_tol*/;
        
        // double tolerance = this->maxValue * Pars::adapt_tol;
        // for (tarLv = baseGrid.maxLevel; tarLv > 0; tarLv--){
        //     if (currValue >= tolerance){
        //         break;
        //     }else{
        //         tolerance *= Pars::adapt_tol;   
        //     }
        // }

        // // Direct Calculation
        // // ******************
        // // Intermediate variable
        // int tarWin = Pars::max_level;     // The target windows
        // // Set up the value of each 
        // for (int d = 0; d < errPredProp.size(); d++){
        //     // [!] Note that every number in propVal is non negative !!!
        //     double propVal = (*errPredProp[d])[_parID];
        //     double normVal = propVal/this->maxValue[d];

        //     // First check (zero bound)
        //     if (normVal < __DBL_EPSILON__){
        //         // tarWin = std::min(tarWin, Pars::max_level);
        //         continue;
        //     }

        //     // Evaluate the current windows level [0-> (max level), 1-> (max_level-1), forth..]
        //     double currWin = std::log(normVal) / this->winTol;
            
        //     // The target windows is the least value
        //     // tarWin = std::min(tarWin, int(currWin));
        //     if (int(currWin) < tarWin){
        //         tarWin = int(currWin);
        //     }
        // }

        // // A prompt to make sure current code is good enough
        // if (tarWin < 0){
        //     ERROR_LOG << "(Log::Adaptation) A target level outside the bound!\n";
        //     throw std::exception();
        // }

        // // Update the target level
        // int tarLv = Pars::max_level - tarWin;

        // Update the TRL of the node container
        // Aliasing to the node where the current particle is contained
        Node *&_node = baseGrid.nodeMap.at(tempPar.nodeID[_parID]);
        // // The node TRL takes the maximum TRL from all particle inside it
        // _node->tarResLv = std::max<int>(_node->tarResLv, tarLv);

        // // // [DEBUG LINE]
        // // TARGET_LVL_LIST[tarLv]++;

        this->targetLevelCalc(_node, _parID, errPredProp);
    }

    // // [DEBUG LINE]
    // for (int _lvl = 0; _lvl < TARGET_LVL_LIST.size(); _lvl++){
    //     std::cout << "The number of particle targeting level " << _lvl << " : " << TARGET_LVL_LIST[_lvl] << "\n";
    // }

    // Flag for adaptation
    #if (N_BODY > 0)
    // #pragma omp parallel for // Cannot be implemented, since particle shares the node ID
    for (size_t _parID = 0; _parID < tempPar.num; _parID++){
        // Set near surface particle to refine at max level
        if (tempPar.isNearSurface[_parID]){
            // Update the TRL of the node container
            // Aliasing to the node where the current particle is contained
            Node *&_node = baseGrid.nodeMap.at(tempPar.nodeID[_parID]);
            // The node containing near surface particle will be refine to the maximum level
            _node->tarResLv = Pars::max_level;
        }
        
        #if (N_BODY != 0 && ADAPT_INSIDE_BODY == 0)
            // Set inside body particle to refine at max level
            if (std::abs(tempPar.chi[_parID] - 1.0) < 0.01){
                // Update the TRL of the node container
                // Aliasing to the node where the current particle is contained
                Node *&_node = baseGrid.nodeMap.at(tempPar.nodeID[_parID]);
                // The node containing near surface particle will be refine to the maximum level
                _node->tarResLv = Pars::max_level;
            }
        #endif
    }
    #endif


    // [!] A custom code section: Add the section to evaluate the active sign
    for (auto &[_ID,_node] : baseGrid.nodeMap){
        if (_node->isLeaf == false) continue;
        // Only set the active mark for leaf node
        bool activeSign = false;
        for (const auto &parID : _node->parList){
            if (tempPar.isActive[parID] == true){
                activeSign = true;
                break;
            }
        }
        _node->isActive = activeSign;
    }

    return;
}

/**
 *  @brief  Set up the target resolution level of each LEAF node.
 *         
 *  @param  baseGrid [UPDATE] The node container which node TRL is to be updated.
 *  @param  tempPar  The template particle to retrieve particle coordinate data.
 *  @param  errPredProp  The container of the selected properties of each particle, sorted by particle ID. (Must be absolute value)
 *  @param  gridCont  The grid container that account all current particle inside the corresponding leaf node.
*/
void GridNodeAdapt::targetLevelv2(GridNode &baseGrid, const Particle &currPar, const std::vector<std::vector<double>*>& errPredProp,
                    std::unordered_map<int,std::vector<int>> &gridCont){
    // In this method, the value of "selected properties" of each particle are being compared to the MAX_VALUE
    // The target resolution level of the given particle takes the matched windows, see the following illustration

    // Set all nodes tarResLv in the GridNode to ROOT prior to the evaluation
    for (auto &[_ID,_node] : baseGrid.nodeMap){
        // Reset the node target level
        _node->tarResLv = ROOT_LEVEL;
    }

    // Section 1
    // *********
    // Evaluate the maximum value of each 
    this->maxValue.resize(errPredProp.size());
    // double max_value [errPredProp.size()];
    for (int d = 0; d < errPredProp.size(); d++){
        this->maxValue[d] = *(std::max_element(errPredProp[d]->begin(), errPredProp[d]->end()));
    }
    
    // Section 2
    // *********
    // Evaluate through all leaf node
    for (auto& [nodeID, parList] : gridCont){
        Node *&_node = baseGrid.nodeMap.at(nodeID);
        // Iterate through all particle
        for (int _parID : parList){
            this->targetLevelCalc(_node, _parID, errPredProp);

            // Additional for near body surface 
            #if (N_BODY > 0)
            if (currPar.isNearSurface[_parID]){
                // Update the TRL of the node container
                // Aliasing to the node where the current particle is contained
                // Node *&_node = baseGrid.nodeMap.at(currPar.nodeID[_parID]);
                // The node containing near surface particle will be refine to the maximum level
                _node->tarResLv = Pars::max_level;
            }
            #endif
        }
    }

    // [!] A custom code section: Add the section to evaluate the active sign
    for (auto &[_ID,_node] : baseGrid.nodeMap){
        if (_node->isLeaf == false) continue;
        // Only set the active mark for leaf node
        bool activeSign = false;
        for (const auto &parID : _node->parList){
            if (currPar.isActive[parID] == true){
                activeSign = true;
                break;
            }
        }
        _node->isActive = activeSign;
    }

    return;
}

/**
 *  @brief  Set up the target resolution level of each LEAF node.
 *         
 *  @param  baseGrid The node container which node TRL is to be updated.
 *  @param  currPar  The template particle to retrieve particle coordinate data.
 *  @param  GridCont  The grid container that account all current particle inside the corresponding leaf node.
*/
void GridNodeAdapt::assignPar2Grid(const GridNode &baseGrid, const Particle &currPar, std::unordered_map<int,std::vector<int>> &GridCont){
    // Section 0
    // *********
    // Assign all particle into the current base Grid
    // AGENT VARIABLE
    // --------------
    // std::vector<std::vector<int>> cntGN;    // Flag to evaluated nodeList
    GridCont.clear();

    // Find the location of the scatter particle
    for (int i = 0; i < currPar.num; i++){
        // Collect common data
        int currLvl = currPar.level[i];
        
        // Define the particle coordinate
        double coord[DIM];
        coord[0] = currPar.x[i];
        coord[1] = currPar.y[i];
        #if (DIM == 3 )
        coord[2] = currPar.z[i];
        #endif

        // [1] FIND THE LOCATION OF THE CURRENT PARTICLE
        int nodeID = baseGrid.pos2ID(coord, currLvl);

        // Check whether current grid is available
        for (int k = currLvl; k > 0; k--){
            // Node at level k -> is the node is existed [?]
            if (baseGrid.nodeMap.count(nodeID) == 0){
                // Node is not existed -> find the parent
                nodeID = baseGrid.findParent(nodeID);
                // _intlSz[i] = rootParSize * std::pow(0.5, k-1);
            }else{
                // Node is existed then the current node will be the leaf node
                break;
            }
        }

        // Check whether current grid is leaf node
        for (int k = currLvl; k < Pars::max_level; k++){
            // Node at level k -> is the node is a leaf node [?]
            if (!baseGrid.nodeMap.at(nodeID)->isLeaf){
                // Node is not a leaf node -> find the corresponding child
                nodeID = baseGrid.pos2ID(coord, k+1);
            }else{
                // Node is existed then the current node will be the leaf node
                break;
            }
        }

        // [2] CREATE THE NODE LIST
        // Put into the node list
        if (GridCont.count(nodeID) == 0){
            // Only put into the queue if the node is not existed
            // nodeList.push_back(nodeID);
            // IDflag.insert({nodeID, true});      // TRUE flag, the node inside "baseGrid"
            // Create the node
            GridCont.insert({nodeID, std::vector<int>{}});
        }
        
        // [3] ASSIGN THE PARTICLE INTO THE NODE LIST
        // Assign the particle into the grid Node container
        GridCont.at(nodeID).push_back(i);
    }
    return;
}

/**
 *  @brief  Target adaptation level (TRL) evaluation of single particle.
 *         
 *  @param  _node [UPDATE] The node which node TRL is to be updated.
 *  @param  parID The particle ID to be evaluated
 *  @param  errPredProp  The container of the selected properties of each particle, sorted by particle ID.
*/
int GridNodeAdapt::targetLevelCalc(Node *_node, int parID, const std::vector<std::vector<double>*>& errPredProp){
    // Direct Calculation
    // ******************
    // Intermediate variable
    int tarWin = Pars::max_level;     // The target windows
    // Set up the value of each 
    for (int d = 0; d < errPredProp.size(); d++){
        // [!] Note that every number in propVal is non negative !!!
        double propVal = (*errPredProp[d])[parID];
        double normVal = propVal/this->maxValue[d];

        // First check (zero bound)
        if (normVal < __DBL_EPSILON__){
            // tarWin = std::min(tarWin, Pars::max_level);
            continue;
        }

        // Evaluate the current windows level [0-> (max level), 1-> (max_level-1), forth..]
        double currWin = std::log(normVal) / this->winTol;
        
        // The target windows is the least value
        // tarWin = std::min(tarWin, int(currWin));
        if (int(currWin) < tarWin){
            tarWin = int(currWin);
        }
    }

    // A prompt to make sure current code is good enough
    if (tarWin < 0){
        ERROR_LOG << "(Log::Adaptation) A target level outside the bound! VAL" << tarWin << "\n";
        throw std::exception();
    }

    // Update the target level
    int tarLv = Pars::max_level - tarWin;

    // Update the TRL of the node container
    // The node TRL takes the maximum TRL from all particle inside it
    _node->tarResLv = std::max<int>(_node->tarResLv, tarLv);

    return tarLv;
}

/**
 *  @brief  Set the flag of node adaptation. Based on the node target adaptation level (TRL),
 *  group the LEAF into three parts "compressList", "refineList", and "idleList".
 *         
 *  @param  baseGrid   The grid node of the evaluation 
*/
void GridNodeAdapt::set_adaptation_flag(GridNode &baseGrid){
    // Iterate through all Node in baseGrid
    for (auto &[_ID,_node] : baseGrid.nodeMap){
        // Only evaluate the LEAF node
        if (_node->isLeaf){
            // Aliasing of the target level and current level of the node
            int &currLv = _node->level;
            int &tarLv = _node->tarResLv;

            // Check whether need adaptation or refinement or not
            if (tarLv < currLv){
                // GROUP IN COMPRESSION LIST
                this->compressList.insert({_ID, true});
                
                // Set up the compression flag
                baseGrid.nodeMap.at(_ID)->needCompression = true;
            }
            else if (tarLv > currLv){
                // GROUP IN REFINEMENT LIST
                this->refineList.insert({_ID, true});
                
                // Duplicate the node into the new grid list
                Node *_dupNode = new Node(_node);
                _dupNode->headNodeID = _dupNode->nodeID;        // Set for adaptation
                this->tempGrid.nodeMap.insert({_ID,_dupNode});
            }
            else{
                // GROUP IN IDLE LIST
                this->idleList.insert({_ID, true});
            }
        }
    }
}

// #pragma endregion

// ===================================================== //
// +--------------- Particle Adaptation ---------------+ //
// ===================================================== //
// #pragma region PARTICLE_UPDATE

/**
 *  @brief  Evaluate the neighbor of the given particle in the list ("parIDList") from
 *  the old particle as the neighbor candidate ("parNghIDList").
 *         
 *  @param  parIDList   List of particle ID which neighbor need to be evaluated.
 *  @param  parNghIDList  List of particle ID as neighbor candidate.
 *  @param  nghSrcPar   The particle source data of neighbor candidate.
*/
void GridNodeAdapt::findParticleNeighbor(const std::vector<int> &parIDList, const std::vector<int> &parNghIDList, const Particle &nghSrcPar){
    // Evaluate all particle in the list
    // printf("Get into the function\n");
    // #pragma omp parallel for         // [IMPORTANT] This paralel make a problem, slow computation
    for (size_t i = 0; i < parIDList.size(); i++){
        // Alias to the current particle to be evaluated
        const auto &evalParID = parIDList[i];
        // printf("Evaluation on particle %d\n", evalParID);

        // Check the distance criteria on all neighbor candidate particle
        for (size_t j = 0; j < parNghIDList.size(); j++){
            // const auto nghParID : parNghIDList
            const auto &nghParID = parNghIDList[j];

            // The internal variable (distance squared)
            double R2 = 0;

            // Calculate the distance square in x direction
            double dx = this->newPar->x[evalParID] - nghSrcPar.x[nghParID];
            dx = dx*dx;
            R2 += dx;
            
            // printf("Evaluation on particle %d (taking x)\n", evalParID);

            // Calculate the distance square in y direction
            if (DIM > 1){
                double dy = this->newPar->y[evalParID] - nghSrcPar.y[nghParID];
                dy = dy*dy;
                R2 += dy;
            }

            // printf("Evaluation on particle %d (taking y)\n", evalParID);

            // Calculate the distance square in z direction
            if (DIM > 2){
                double dz = this->newPar->z[evalParID] - nghSrcPar.z[nghParID];
                dz = dz*dz;
                R2 += dz;
            }

            // printf("Evaluation on particle %d (taking z)\n", evalParID);

            // Check neighbor distance criteria (Using the base size of old data)
            // Support radius ('s' is the size of Head particle)
            double supRad = this->newPar->s[evalParID] * Pars::r_sup * Pars::r_buff;
            
            // Remember checking the squared radius
            if (R2 < supRad*supRad){
                // Inside the support radius, put this particle inside the particle neighbor
                this->newPar->neighbor[evalParID].push_back(nghParID);

                // if (evalParID == 0){
                //     std::cout << "GOOD DAY BUDDY\n";
                //     NUM_STR++;
                // }
            }
        }
        
        // printf("Done Evaluation on particle %d\n", evalParID);
    }

    return;
}

/**
 *  @brief  Generate the particle in the given node. 
 *  NOTE: The size of particle is using the size of particle in the head node.
 *         
 *  @param  currNode  The node where the particle are to be generated.
 *  @param  headSize  The size of particle in the head node.
*/
void GridNodeAdapt::generateParticle(Node *currNode, double _headSize){
    // Generate the particle in the current Node

    // Update the current node data
    currNode->parList.clear();  // Reserve the particle ID container on the current node 
    currNode->headNodeID = -1;  // Turn off the head node ID flag

    // Size of the particle inside the node
    double _parSize = currNode->length / this->parNum;
    // double _headSize = headNode->length / this->parNum;

    // Generate all particle inside the current Node
    int locParIndex[DIM];     // The particle local index (taken from the node)
    for (int locParID = 0; locParID < this->totParNum; locParID++){
        // Put the current particle ID into the node
        currNode->parList.push_back(this->newPar->num);    // Push the last pushed ID inside the parList

        // The local index coordinate inside the node
        basis_loop(d) locParIndex[d] = (locParID/this->div[d]) % this->parNum;

        // Assign the particle position
        // The x position
        double _x = currNode->pivCoor[0] + (0.5 + locParIndex[0])*_parSize;
        this->newPar->x.push_back(_x);
        
        // The y position
        double _y = currNode->pivCoor[1] + (0.5 + locParIndex[1])*_parSize;
        this->newPar->y.push_back(_y);
        
        // The z position
        if (DIM > 2) {
            double _z = currNode->pivCoor[2] + (0.5 + locParIndex[2])*_parSize;
            this->newPar->z.push_back(_z);
        }

        // Assign other data
        this->newPar->nodeID.push_back(currNode->nodeID);  // Must be targeting the real node ID
        this->newPar->level.push_back(currNode->level);    // Also the real level
        this->newPar->s.push_back(_headSize);    // Size using the old size
        this->newPar->isActive.push_back(currNode->isActive);    // Flag will be set as true (because if divided it means an active block)
        this->newPar->neighbor.push_back({});    // Empty neighbor list
        this->newPar->num++;     // Add the number of the particle
    }
    return;
}

// /**
//  *  @brief  Generate the particle in the given node. 
//  *  NOTE: The size of particle is using the size of particle in the head node.
//  *         
//  *  @param  currNode  The node where the particle are to be generated.
//  *  @param  headSize  The size of particle in the head node.
//  *  @param  activeSign  The active sign for the generated particle.
// */
// void GridNodeAdapt::generateParticleActive(Node *currNode, double _headSize, bool _isActive){
//     // Generate the particle in the current Node

//     // Update the current node data
//     currNode->parList.clear();  // Reserve the particle ID container on the current node 
//     currNode->headNodeID = -1;  // Turn off the head node ID flag

//     // Size of the particle inside the node
//     double _parSize = currNode->length / this->parNum;
//     // double _headSize = headNode->length / this->parNum;

//     // Generate all particle inside the current Node
//     int locParIndex[DIM];     // The particle local index (taken from the node)
//     for (int locParID = 0; locParID < this->totParNum; locParID++){
//         // Put the current particle ID into the node
//         currNode->parList.push_back(this->newPar->num);    // Push the last pushed ID inside the parList

//         // The local index coordinate inside the node
//         basis_loop(d) locParIndex[d] = (locParID/this->div[d]) % this->parNum;

//         // Assign the particle position
//         // The x position
//         double _x = currNode->pivCoor[0] + (0.5 + locParIndex[0])*_parSize;
//         this->newPar->x.push_back(_x);
        
//         // The y position
//         double _y = currNode->pivCoor[1] + (0.5 + locParIndex[1])*_parSize;
//         this->newPar->y.push_back(_y);
        
//         // The z position
//         if (DIM > 2) {
//             double _z = currNode->pivCoor[2] + (0.5 + locParIndex[2])*_parSize;
//             this->newPar->z.push_back(_z);
//         }

//         // Assign other data
//         this->newPar->nodeID.push_back(currNode->nodeID);   // Must be targeting the real node ID
//         this->newPar->level.push_back(currNode->level);     // Also the real level
//         this->newPar->s.push_back(_headSize);           // Size using the old size
//         this->newPar->isActive.push_back(_isActive);    // Flag will be set as true (because if divided it means an active block)
//         this->newPar->neighbor.push_back({});           // Empty neighbor list
//         this->newPar->num++;     // Add the number of the particle
//     }
//     return;
// }

// #pragma endregion


/**
 *  @brief  Perform a particle list division from the current node to it's child-LEAF node.
 *
 *  @param  currNode     Node which the particle list to be divided.
 *  @param  baseParticle Basis particle object to provide individual particle coordinate
 *  for the calculation.
 *  @param  chdIDlist    List of all child IDs.
 *  @param  type    Type of particle division, [1] REFINEMENT, [DEFAULT] COMPRESSION
 * 
 *  @return A check mark of 0 or 1 whether the node is refined (1) or not (0).
*/
int GridNodeAdapt::assignParticle(Node *currNode, const Particle &templatePar, const std::vector<int> &chdIDlist, int type){
    // // Give a check, check whether all child Node is valid
    // // bool invalid = true;
    // for (auto &chdID : chdIDlist){
    //     // No target child
    //     if (this->nodeMap.count(chdID) == 0){
    //         ERROR_LOG << "Node " << chdID << " as child of node " << currNode->nodeID << " is not existed!\n";
    //         throw std::runtime_error("ERROR [GRID NODE] Unable to get the child, the node is not existed in the grid!");
    //     }

    //     // Reserve the child particle list data
    //     this->nodeMap.at(chdID)->parList.clear();
        
    //     // Additional for parameter update for adaptation (REFINEMENT)
    //     if (type == 1) this->nodeMap.at(chdID)->headNodeID = currNode->headNodeID;

    //     // // Child list is invalid 
    //     // if (chdID != 0){
    //     //     invalid = false;
    //     // }
    // }

    // // // Error Handling
    // // if (invalid){
    // //     throw std::runtime_error("ERROR [GRID NODE] Invalid child ID list!");
    // // }

    // // Internal variable
    // double chdNodeLength = currNode->length / 2.0;
    // int posIdx[DIM];        // Temporary child node index container
    // double parPos[DIM];     // Temporary particle coordinate holder
    
    // // Move the particle into the child node
    // for (int &_particleID : currNode->parList){
    //     // Retrieve the position of the particle
    //     parPos[0] = templatePar.x[_particleID];
    //     if (DIM > 1) parPos[1] = templatePar.y[_particleID];
    //     if (DIM > 2) parPos[2] = templatePar.z[_particleID];
    //     basis_loop(d){
    //         posIdx[d] = (parPos[d] - currNode->pivCoor[d]) / chdNodeLength;
    //         // Undefined behaviour handler
    //         if (posIdx[d] < 0 || posIdx[d] > 1){
    //             // The local index is out of range
    //             ERROR_LOG << "Particle " << _particleID << " is outside the current node " << currNode->nodeID << " box! (" << posIdx[d] << "," << (parPos[d] - currNode->pivCoor[d]) << "," << chdNodeLength<< ")\n";
    //             throw std::range_error("ERROR [GRID NODE] Child index is out of range!");
    //         }
    //     }
        
    //     // Find the child local index
    //     int chdLocID = 0;
    //     // This following line is where the previous wrong algorithm occur [The line issue has resolved]
    //     basis_loop(d) chdLocID = chdLocID*2 + posIdx[(DIM-1) - d];

    //     // Copy the particle ID into the corresponding child Node
    //     Node *&_chdNode = this->nodeMap.at(chdIDlist[chdLocID]);
    //     _chdNode->parList.push_back(_particleID);
    // }

    // // Additional for parameter update for adaptation (REFINEMENT)
    // if (type == 1){
    //     // Release the adaptation (REFINMENT) parameter of the current node
    //     currNode->headNodeID = -1;
    //     currNode->parList.clear();
    // }

    return 1;
}