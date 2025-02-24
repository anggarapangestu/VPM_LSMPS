#include "neighbor.hpp"

/**
 *  @brief A grid generator, used as a temporary data grouping for neighbor evaluation.
 *  NOTE: Works well both for 2D and 3D simulation.
 *  
 *  @param	_parEval  The particle data for grid construction.
 *  @param	_tempGrid  [OUTPUT] The temporary grid to be generated.
 */
void neighbor::create_temp_grid(const Particle& _par, NghBaseGrid &_grid){
    // Evaluate extreme coordinate
	double maxCoor[DIM];
	double minCoor[DIM];
	// Evaluate extreme in x direction
		minCoor[0] = *(std::min_element(_par.x.begin(),_par.x.end()));
		maxCoor[0] = *(std::max_element(_par.x.begin(),_par.x.end()));
	// Evaluate extreme in y direction
		minCoor[1] = *(std::min_element(_par.y.begin(),_par.y.end()));
		maxCoor[1] = *(std::max_element(_par.y.begin(),_par.y.end()));
	// Evaluate extreme in z direction
	if (DIM > 2){
		minCoor[2] = *(std::min_element(_par.z.begin(),_par.z.end()));
		maxCoor[2] = *(std::max_element(_par.z.begin(),_par.z.end()));
	}
    // Perform grid initialization
    _grid.initialize_grid(maxCoor, minCoor, _par.s);

    return;
}

/**
 *  @brief Neighbor search function manager for neighbor search toward 
 *  the particle data itself.
 *  NOTE: Works well both for 2D and 3D simulation.
 *  
 *  @param	_evalPar  The particle data for neighbor evaluation.
 *  @param  _baseGrid  The grid node data for Grid Node type neighbor search only.
 */
void neighbor::neigbor_search(Particle &parEval, const GridNode &_baseGrid){
    // Print log to console
    printf("Evaluating neighbor ...\n");
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif


    // Direct neighbor search evaluation
    if (Pars::opt_neighbor == 0){
        printf("%s<+> Type 0: Direct neighbor search %s\n", FONT_CYAN, FONT_RESET);
        
        // Evaluate the neighbor
        directFindNgh _directFind;
        _directFind.find_neighbor(parEval.neighbor, parEval.s, 
                                 parEval.x, parEval.y, parEval.z);
    }
    
    // Link list neighbor search evaluation
    else if (Pars::opt_neighbor == 1){
        printf("%s<+> Type 1: Link list neighbor search %s\n", FONT_CYAN, FONT_RESET);
        
        // Create temporary grid
        NghBaseGrid _tempGrid;
        this->create_temp_grid(parEval, _tempGrid);

        // Evaluate the neighbor
        LinkListNgh _linkList;
        _linkList.find_neighbor(parEval.neighbor, _tempGrid, parEval.s,
                                parEval.x, parEval.y, parEval.z);
    }
    
    // [NEED FURTHER REFACTOR]
    // Included neighbor search package from Cell List
    else if (Pars::opt_neighbor == 2){
        printf("%s<+> Type 2: Cell list neighbor search %s\n", FONT_CYAN, FONT_RESET);
        this->cellListData.findNeighbor(parEval);
        // this->cellListData.checkNGH(parEval);
    }

    // Spatial hash neighbor search evaluation
    else if (Pars::opt_neighbor == 3){
        printf("%s<+> Type 3: Spatial hash neighbor search %s\n", FONT_CYAN, FONT_RESET);
        
        // Create temporary grid
        NghBaseGrid _tempGrid;
        this->create_temp_grid(parEval, _tempGrid);

        // Evaluate the neighbor
        SpatialHashNgh _spatialHash;
        _spatialHash.find_neighbor(parEval.neighbor, _tempGrid, parEval.s,
                                   parEval.x, parEval.y, parEval.z);
    }

    // Grid node neighbor search evaluation
    else if (Pars::opt_neighbor == 4){
        printf("%s<+> Type 4: Grid node neighbor search %s\n", FONT_CYAN, FONT_RESET);
        
        // Evaluate the neighbor
        GridNodeNgh _gridNode;
        _gridNode.find_neighbor(parEval.neighbor, _baseGrid, parEval);
    }
    
    // Neighbor search summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Neighbor search computation time:  [%f s]\n", _time);
}

/**
 *  @brief Neighbor search function manager for two different particle data set.
 *  Evaluation is done in one direction to evaluate the source particle data as
 *  the neighbor of the target particle data.
 *  NOTE: Works well both for 2D and 3D simulation.
 *  
 *  @param	_evalPar  The particle target particle data for neighbor evaluation.
 *  @param	_srcPar  The source particle data, work as a data format.
 */
void neighbor::neigbor_search(Particle& _evalPar, Particle &_srcPar, std::vector<std::vector<int>> &_nghIDList){
    // Print log to console
    printf("Evaluating inter search neighbor ...\n");
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // Evaluate the neighbor
    InterSearchNgh _interSearch;
    _interSearch.find_neighbor(_srcPar, _evalPar, _nghIDList);
    
    // Neighbor search summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Neighbor search computation time:  [%f s]\n", _time);
}


// ************************************************************************************
// ====================================================================================

// [OLD PACKAGE] The particle neighbor list
void neighbor::cell_list_init(Particle& parEval){
    // Generate the Cell List
    printf("Generating cell list ...\n");
    clock_t t = clock();

    // Initialize the Cell List
    cellListData.initCellList(parEval);
    
    // Evaluate the particle cell ID
    cellListData.createCellList(parEval);

    // // Show particle neighbor
    // cellListData.showBasisNgh();
    
    // particle generation initialization summary time display
    t = clock() - t;
	printf("<-> Cell list initialization\n");
    printf("    comp. time:                        [%f s]\n", (double)t/CLOCKS_PER_SEC);

    // Evaluate the particle cell ID
    if (Pars::flag_save_cell){cellListData.saveCellData();}
    
}

// ************************************************************************************
// ====================================================================================

// [OLD PACKAGE] Particle adaptation 
bool neighbor::particle_adaptation(const Particle& parEval, Particle& baseParticle, std::vector<double>& PARsize){
    // Note: At initial the distribution of parEval and baseParticle is the same but not the properties

    // computational time accumulation
    clock_t t;
    t = clock();
    /* Procedure of particle adaptation:
       [1] Evaluate each particle, determine if need adaptation of not
       [2] Perform the particle adaptation if necesarry (based on the procedure 1)
    */

    // PROCEDURE 1 ! : Evaluating particle adaptation
    // *************
    // The marker whether the adaptation need to performed or not
    bool _adaptation = false;

    // TODO: Find maximum vorticity
    double vor_max = 0.0e0;
    for (int i = 0; i < parEval.num; i++)
    {
        vor_max = vor_max > std::abs(parEval.vorticity[i]) ? 
                  vor_max : std::abs(parEval.vorticity[i]);
    }

    // Initialize the Cell Level Up setting (Setting for particle adaptation)
    this->cellListData.setLevelUp(0, 0, 0);

    // Assign the cell for Level Up (evaluate each particle)
    for (int i = 0; i < parEval.num; i++)
    {   
        // ... LEVEL 1 CHECK ...
        // Adaptation evaluation still based on vorticity value of the evaluated particle
        if (std::abs(parEval.vorticity[i]) >= 1.0e-5 * vor_max)  // Threshold set to be 1.0e-5 of the max value,
                                                                 //   set lower than the particle redistribution
        {
            // ... LEVEL 2 CHECK ...
            // Do adaptation if the particle level != maxlevel (or lower) and the
            if (parEval.level[i] < Pars::max_level){
                // Assign the cell for Level Up
                this->cellListData.setLevelUp(parEval.basis_label[i], parEval.cell_label[i], 1);

                // Do the adaptation
                _adaptation = true;
            }
        }
    }
    
    // Adaptation procedure:
    // -> List all the cell to be devided (save the basis-cell-ID pair)
    // -> Find the neighboring cell to be divided (save the basis-cell-ID pair)
    // -> At each cell: divide the particle into the target level (only if par level < cell level)
    // -> Put the particle into the corresponding cell
    // -> Devide the cell
    
    // PROCEDURE 2 ! : Performing particle adaptation
    // *************
    if (_adaptation == true)
    {
        // Perform the adaptive algorithm
        this->cellListData.performAdaptive(baseParticle, 0, Pars::max_level, PARsize);

        // Particle adaptation summary time display
        t = clock() - t;
        printf("<+> Particle adaptation done ...\n");
        printf("<+> Number of particle after adaptation : %8d\n", baseParticle.num);
        printf("<-> Particle adaptation calculation \n");
        printf("    comp. time:                        [%f s]\n\n", (double)t/CLOCKS_PER_SEC);
        
        // Evaluate the verlet list of new particle distribution
        // this->neigbor_search(baseParticle);		// Already taken on the particle adaptation
    }else{
        printf("<+> No particle adaptation\n");
    }

    return _adaptation;
}

/**
 *  @brief Neighbor search function manager for two different particle data set.
 *  Evaluation is done in one direction whether to find neighbor of scatter particle
 *  containing grid particle or other way around. The calculation direction is dictated 
 *  by the input direction. The nodegrid containing the grid particle.
 *  NOTE: Works well both for 2D and 3D simulation.
 *  
 *  @param	_gridPar  The particle lies on the grid, arranged by the grid Node.
 *  @param	_scatPar  The scatter particle data.
 *  @param	_baseGrid  The grid node arranging the grid particle.
 *  @param	_dir  The neighbor calculation direction.
 *                   > 1:= Neighbor of scatter by using grid;
 *                   > 2:= Neighbor of grid by using scatter
 *  @param	_nghIDList  The storage for neighbor ID list [OUTPUT -> Used for LSMPS]
 *                   :: The indexing [currParID][nghParID]
 *                   >  _dir = 1 --> _nghIDList is the SCATTER particle neighbor
 *                   >  _dir = 2 --> _nghIDList is the GRID particle neighbor
 *  @param	_intlSz  The size of updated particle (OUTPUT used for LSMPS)
 */
void neighbor::inter_neigbor_search(Particle &_gridPar, 
		Particle &_scatPar, const GridNode &_baseGrid, 
		int _dir, std::vector<std::vector<int>> &_nghIDList, 
		std::vector<double> &_intlSz){
    // Procedure
    // [1] Create a new container of grid block containing the scatter particle. (Based on the direction)
    // [2] Inter neighbor calculation.

    MESSAGE_LOG << "CALCULATING THE GRID NODE INTER-NEIGHBOR!\n";

    // Variable aliasing 
    // *****************
    Particle& _gP = _gridPar;
    Particle& _sP = _scatPar;
    const GridNode& _nG = _baseGrid;
    double rootParSize = Pars::sigma * std::pow(2.0, Pars::max_level);

    // Define the neighbor list
    // std::vector<std::vector<int>> *_nghIDPtr;       // A dummy pointer container
    
    // INITIALIZATION 
    // **************
    // Reserve the data container (for output)
    _nghIDList.clear();
    _intlSz.clear();  // The size of the scattered particle
    if (_dir == 1){     
        // // Evaluating the NGH of SCATTER particle
        // _scatPar.neighbor.clear();
        // _scatPar.neighbor.resize(_scatPar.num);
        // _nghIDPtr = &_scatPar.neighbor;
        
        // Update the size of particle
        _nghIDList.resize(_sP.num);
        _intlSz.resize(_sP.num);
    }else if (_dir == 2){
        // // Evaluating the NGH of GRID particle
        // _gridPar.neighbor.clear();
        // _gridPar.neighbor.resize(_gridPar.num);
        // _nghIDPtr = &_gridPar.neighbor;

        // Update the size of particle
        _nghIDList.resize(_gP.num);
        _intlSz.resize(_gP.num);
    }
    // std::vector<std::vector<int>> &_nghIDList = *_nghIDPtr;

    
    // NEIGHBOR EVALUATION
    // *******************
    if(_dir == 1){
    // ===================================================
    // ------------------- DIRECTION 1 -------------------
    // ===================================================
    // Direction 1 -> Calculating the scatter particle neighbor

    // PROCEDURE:
    //  1. Assign the SCATTER particle into the grid node at very leaf node
    //     - Sub 1: At each SCATTER particle find the leaf node on the location
    //     - Sub 2: Create a node container that contain all SCATTER particle in current node
    //     - Sub 3: Assign the interpolation particle size corresponding to node size
    //  2. Find neighbor of SCATTER particle toward the GRID particle
    //     - Sub 1: Find all leaf node as a neighbor coresponding to current node
    //     - Sub 2: Evaluate the neighbor on the grid particle located on ngh node

    // AGENT VARIABLE
    // --------------
    // std::vector<std::vector<int>> cntGN;    // Flag to evaluated nodeList
    std::unordered_map<int,std::vector<int>> cntGN;    // Flag to evaluated nodeList
    std::vector<int> nodeList;              // List of all node containing scattered element (Actually the collection of all leaf node)
    // std::unordered_map<int,int> _map;    // Flag to evaluated nodeList
    // WARNING_LOG << "Assigning the scatter particle container\n";

    // PROCEDURE 1
    // ```````````
    // Find the location of the scatter particle
    for (int i = 0; i < _sP.num; i++){
        // Collect common data
        int currLvl = _sP.level[i];

        // Update the interpolation size (initially similar toward the current particle)
        // _intlSz[i] = _sP.s[i];
        _intlSz[i] = rootParSize * std::pow(0.5, currLvl);
        
        // Define the particle coordinate
        double coord[DIM];
        coord[0] = _sP.x[i];
        coord[1] = _sP.y[i];
        #if (DIM == 3 )
        coord[2] = _sP.z[i];
        #endif

        // [1] FIND THE LOCATION OF THE CURRENT PARTICLE
        int nodeID = _nG.pos2ID(coord, currLvl);

        // Check whether current grid is available
        for (int k = currLvl; k > 0; k--){
            // Node at level k -> is the node is existed [?]
            if (_nG.nodeMap.count(nodeID) == 0){
                // Node is not existed -> find the parent
                nodeID = _nG.findParent(nodeID);
                _intlSz[i] = rootParSize * std::pow(0.5, k-1);
            }else{
                // Node is existed then the current node will be the leaf node
                // if (!_nG.nodeMap.at(nodeID)->isLeaf){
                //     throw std::runtime_error("ERROR [INTER NEIGHBOR] : First available parent is not a leaf!");
                // }
                break;
            }
        }

        // Check whether current grid is leaf node
        for (int k = currLvl; k < Pars::max_level; k++){
            // Node at level k -> is the node is a leaf node [?]
            if (!_nG.nodeMap.at(nodeID)->isLeaf){
                // Node is not a leaf node -> find the corresponding child
                nodeID = _nG.pos2ID(coord, k+1);
                _intlSz[i] = rootParSize * std::pow(0.5, k+1);
            }else{
                // Node is existed then the current node will be the leaf node
                break;
            }
        }

        // [2] CREATE THE NODE LIST
        // Put into the node list
        if (cntGN.count(nodeID) == 0){
            // Only put into the queue if the node is not existed
            nodeList.push_back(nodeID);
            // IDflag.insert({nodeID, true});      // TRUE flag, the node inside "baseGrid"
            // Create the node
            cntGN.insert({nodeID, std::vector<int>{}});
        }
        
        // [3] ASSIGN THE PARTICLE INTO THE NODE LIST
        // Assign the particle into the grid Node container
        cntGN.at(nodeID).push_back(i);
    }

    
    // PROCEDURE 2
    // ```````````
    // Evaluating Neighbor
	// Iterate through all SCATTER Particle in the grid
	#pragma omp parallel for			// [-> Create a problem in the server computer [RESOLVED by create a container 'leafNodeList' first]]
	for (size_t n = 0; n < nodeList.size(); n++){
	// for (const auto&[_nodeID,_currNode] : _baseGrid.nodeMap){
		// Alias to the current node and particle list
        int _nodeID = nodeList.at(n);
		const Node* _currNode = _nG.nodeMap.at(_nodeID);
        std::vector<int>& currParList = cntGN.at(_nodeID);

        // MESSAGE_LOG << "Node " << _currNode->nodeID << "\n";
        
        // Get the list of neighbor node
        std::vector<Node*> _nghNodeList;
        _nG.findNghAll(_nghNodeList, _currNode);

        // std::cout << "> Done node ngh finding of " << _currNode->nodeID << "\n";

		// Iterate through all scatter particle inside the current node (it will share the same neighbor node)
		for (size_t i = 0; i < currParList.size(); i++){
            // Aliasing the current particle ID (scattered particle)
			const int &ID_i = currParList.at(i);

            // std::cout << ">> Particle " << ID_i << "\n";
			
            // Evaluate to each point inside the grid neighbor
			for (const auto &nghNode : _nghNodeList){
                // if (cntGN.count(nghNode->nodeID) == 0) continue;
				for (const auto &ID_j : nghNode->parList){		
					// // An exception for include itself on neighbor list
					// if (!Pars::flag_ngh_include_self && (ID_j == ID_i)) continue;

                    // Calculate the distance square between two points
					double _dr2 = 0;
					// Calculate in x direction
						double _dx = _sP.x[ID_i] - _gP.x[ID_j];
						_dr2 += (_dx*_dx);
					// Calculate in y direction
						double _dy = _sP.y[ID_i] - _gP.y[ID_j];
						_dr2 += (_dy*_dy);
					// Calculate in z direction
					if (DIM > 2){
						double _dz = _sP.z[ID_i] - _gP.z[ID_j];
						_dr2 += (_dz*_dz);
					}

					// Evaluate distance
					if (Pars::opt_ngh_interact == 1){
						// Check the particle i
						double _rSup = _intlSz[ID_i] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
						if (_dr2 < (_rSup*_rSup)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
					else if (Pars::opt_ngh_interact == 2){
						// Support radius of average size
						double _rSupAve = ((_intlSz[ID_i] + _gP.s[ID_j]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
						if (_dr2 < (_rSupAve*_rSupAve)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
				}
			}
		}
        // std::cout << "> Done ngh evaluation of " << _currNode->nodeID << "\n";
	}
    }
    else if(_dir == 2){
    // ===================================================
    // ------------------- DIRECTION 2 -------------------
    // ===================================================
    // Direction 2 -> Calculating the neighbor of grid particle relative to scattered particle

    // PROCEDURE:
    //  1. Assign the SCATTER particle into the grid node at very leaf node
    //     - Sub 1: At each SCATTER particle find the leaf node on the location
    //     - Sub 2: Create a node container that contain all SCATTER particle in current node
    //     - Sub 3: Assign the interpolation particle size corresponding to node size -> Grid node
    //              - Set the evaluation level of the current node (Not the orginal level)
    //  2. Find neighbor of SCATTER particle toward the GRID particle
    //     - Sub 1: Find all leaf node as a neighbor coresponding to current node
    //     - Sub 2: Evaluate the neighbor on the grid particle located on ngh node

    // AGENT VARIABLE
    // --------------
    // std::vector<std::vector<int>> cntGN;    // Flag to evaluated nodeList
    std::unordered_map<int,std::vector<int>> cntGN;    // Flag to evaluated nodeList (containing the scatter particle)
    std::unordered_map<int,int> tarLvl;        // The container of current node ngh evaluation level
    std::vector<int> nodeList;         // List of all node containing scattered element
    // std::vector<int> tarLvl;           // List of all node containing scattered element

    // Find the location of the scatter particle
    for (int i = 0; i < _sP.num; i++){
        // Collect common data
        int currLvl = _sP.level[i];
        
        // Define the particle coordinate
        double coord[DIM];
        coord[0] = _sP.x[i];
        coord[1] = _sP.y[i];
        #if (DIM == 3 )
        coord[2] = _sP.z[i];
        #endif

        // [1] FIND THE LOCATION OF THE CURRENT PARTICLE
        int nodeID = _nG.pos2ID(coord, currLvl);

        // Check whether current grid is available
        for (int k = currLvl; k > 0; k--){
            // Node at level k -> is the node is existed [?]
            if (_nG.nodeMap.count(nodeID) == 0){
                // Node is not existed -> find the parent
                nodeID = _nG.findParent(nodeID);
                // _intlSz[i] = rootParSize * std::pow(0.5, k-1);
            }else{
                // Node is existed then the current node will be the leaf node
                // if (!_nG.nodeMap.at(nodeID)->isLeaf){
                //     throw std::runtime_error("ERROR [INTER NEIGHBOR] : First available parent is not a leaf!");
                // }
                break;
            }
        }

        // Check whether current grid is leaf node
        for (int k = currLvl; k < Pars::max_level; k++){
            // Node at level k -> is the node is a leaf node [?]
            if (!_nG.nodeMap.at(nodeID)->isLeaf){
                // Node is not a leaf node -> find the corresponding child
                nodeID = _nG.pos2ID(coord, k+1);
                // _intlSz[i] = rootParSize * std::pow(0.5, k+1);
            }else{
                // Node is existed then the current node will be the leaf node
                break;
            }
        }

        // [2] CREATE THE NODE LIST
        // Put into the node list
        if (cntGN.count(nodeID) == 0){
            // Only put into the queue if the node is not existed
            nodeList.push_back(nodeID);
            // IDflag.insert({nodeID, true});      // TRUE flag, the node inside "baseGrid"
            // Create the node
            cntGN.insert({nodeID, std::vector<int>{}});
            tarLvl.insert({nodeID, currLvl});
        }
        
        // [3] ASSIGN THE PARTICLE INTO THE NODE LIST
        // Assign the particle into the grid Node container
        cntGN.at(nodeID).push_back(i);
        tarLvl[nodeID] = std::min(tarLvl[nodeID], currLvl); // Set the level of largest particle (which is the minimum level)
    }

    // Evaluating Neighbor
    // *******************
    // Collect all leaf node only (for the sake of parallel programming)
	std::vector<const Node*> leafNodeList;
	for (const auto&[_nodeID,_currNode] : _nG.nodeMap){
		if(!(_currNode->isLeaf)) continue;
		leafNodeList.push_back(_currNode);
	}
	
	// Iterate through all node in the grid
	// #pragma omp parallel for			// [-> Create a problem in the server computer [RESOLVED by create a container 'leafNodeList' first]]
	for (size_t n = 0; n < leafNodeList.size(); n++){
	// for (const auto&[_nodeID,_currNode] : _baseGrid.nodeMap){
		// Alias to the current node
        // MESSAGE_LOG << "Start Iteration\n";
		const Node *&_currNode = leafNodeList[n];
        // std::cout << "Done taking the node\n";
        int nodeID = _currNode->nodeID;
        // std::cout << "Done taking the ID: " << nodeID << "\n";
        int nghLvTar = tarLvl.at(nodeID);
        // std::cout << "Done taking thelevel\n";
        // std::cout << "Evaluating the node " << nodeID << "\n";
        
        // Get the list of neighbor node
        std::vector<Node*> _nghNodeList;
        _nG.findNghAll(_nghNodeList, _currNode, nghLvTar);

        // std::cout << "> Done node ngh finding of " << _currNode->nodeID << "\n";

		// Iterate through all point inside the current node (it will share the same neighbor node)
		for (size_t i = 0; i < _currNode->parList.size(); i++){
			// Aliasing the current particle ID
			const int &ID_i = _currNode->parList[i];
            _intlSz[ID_i] = rootParSize * std::pow(0.5, nghLvTar);

            // std::cout << ">> Particle " << ID_i << "\n";
			
            // Evaluate to each point inside the grid neighbor
			for (const auto &nghNode : _nghNodeList){
                if (cntGN.count(nghNode->nodeID) == 0) continue;
				for (const auto &ID_j : cntGN.at(nghNode->nodeID)){		
					// // An exception for include itself on neighbor list
					// if (!Pars::flag_ngh_include_self && (ID_j == ID_i)) continue;

                    // Calculate the distance square between two points
					double _dr2 = 0;
					// Calculate in x direction
						double _dx = _gP.x[ID_i] - _sP.x[ID_j];
						_dr2 += (_dx*_dx);
					// Calculate in y direction
						double _dy = _gP.y[ID_i] - _sP.y[ID_j];
						_dr2 += (_dy*_dy);
					// Calculate in z direction
					if (DIM > 2){
						double _dz = _gP.z[ID_i] - _sP.z[ID_j];
						_dr2 += (_dz*_dz);
					}

					// Evaluate distance
					if (Pars::opt_ngh_interact == 1){
						// Check the particle i
						double _rSup = _intlSz[ID_i] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
						if (_dr2 < (_rSup*_rSup)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
					else if (Pars::opt_ngh_interact == 2){
						// Support radius of average size
						double _rSupAve = ((_intlSz[ID_i] + _sP.s[ID_j]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
						if (_dr2 < (_rSupAve*_rSupAve)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
				}
			}
		}
        // std::cout << "> Done ngh evaluation of " << nodeID << "\n\n";
	}
    }

    return;
}