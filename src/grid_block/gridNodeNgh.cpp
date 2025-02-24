#include "gridNodeNgh.hpp"

/**
 *  @brief  Evaluate the neighboring interaction for each particle. 
 *  Interaction pairs are determined by grid node relation.
 *
 *  @param  _nghIDList	[OUTPUT] Neighbor ID list.
 *  @param  _baseGrid 	The node container as neighbor evaluation tools.
 *  @param  _evalPar  	The particle to evaluate neighbor.
*/
void GridNodeNgh::find_neighbor(std::vector<std::vector<int>> &_nghIDList, 
								const GridNode &_baseGrid, const Particle &_par)
{
    // **Evaluate neighbor interation
    _nghIDList.clear();
    _nghIDList.resize(_par.num);
    std::vector<bool> evalFlag(_par.num, false);

	// Collect all leaf node only (for the sake of parallel programming)
	std::vector<const Node*> leafNodeList;
	for (const auto&[_nodeID,_currNode] : _baseGrid.nodeMap){
		if(!(_currNode->isLeaf)) continue;
		leafNodeList.push_back(_currNode);
	}
	
	// Iterate through all node in the grid
	#pragma omp parallel for			// [-> Create a problem in the server computer [RESOLVED by create a container 'leafNodeList' first]]
	for (size_t n = 0; n < leafNodeList.size(); n++){
	// for (const auto&[_nodeID,_currNode] : _baseGrid.nodeMap){
		// Alias to the current node
		const Node *&_currNode = leafNodeList[n];

		// // Only evaluate the leaf node
        // if(!(_currNode->isLeaf)) continue;
        
        // Get the list of neighbor node
        std::vector<Node*> _nghNodeList;
        _baseGrid.findNghAll(_nghNodeList, _currNode);

		// Iterate through all point inside the current node (it will share the same neighbor node)
		for (size_t i = 0; i < _currNode->parList.size(); i++){
			// Aliasing the current particle ID
			const int &ID_i = _currNode->parList[i];
			
            // Evaluate to each point inside the grid neighbor
			for (const auto &nghNode : _nghNodeList){
				for (const auto &ID_j : nghNode->parList){		
					// An exception for include itself on neighbor list
					if (!Pars::flag_ngh_include_self && (ID_j == ID_i)) continue;

                    // Calculate the distance square between two points
					double _dr2 = 0;
					// Calculate in x direction
						double _dx = _par.x[ID_i] - _par.x[ID_j];
						_dr2 += (_dx*_dx);
					// Calculate in y direction
						double _dy = _par.y[ID_i] - _par.y[ID_j];
						_dr2 += (_dy*_dy);
					// Calculate in z direction
					if (DIM > 2){
						double _dz = _par.z[ID_i] - _par.z[ID_j];
						_dr2 += (_dz*_dz);
					}

					// Evaluate distance
					if (Pars::opt_ngh_interact == 1){
						// Check the particle i
						double _rSup = _par.s[ID_i] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
						if (_dr2 < (_rSup*_rSup)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
					else if (Pars::opt_ngh_interact == 2){
						// Support radius of average size
						double _rSupAve = ((_par.s[ID_i] + _par.s[ID_j]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
						if (_dr2 < (_rSupAve*_rSupAve)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
				}
			}
		}
	}
    return;
}


/**
 *  @brief  Evaluate the neighboring interaction for each particle. 
 *  Interaction pairs are determined by grid node relation.
 *
 *  @param  _nghIDList	[OUTPUT] Neighbor ID list of eval particle toward source particle.
 *  @param  _evalPar 	The particle to evaluate neighbor. (The Target Particle)
 *  @param  _evalParNodeMap The evaluated particle mapping.
 *  @param  _baseGrid  	The node container of source particle.
 *  @param  _sourcePar	The source particle to evaluate neighbor. (Properties particle)
*/
void GridNodeNgh::find_inter_neighbor(
	std::vector<std::vector<int>> &_nghIDList,
	const Particle &_par,
	const std::unordered_map<int, std::vector<int>> &_parMap,
	const GridNode &_baseGrd,
	const Particle &_srcPar)
{
	// **Evaluate neighbor interation (spatial hash)
    _nghIDList.clear();
    _nghIDList.resize(_par.num);

	// Collect all leaf node only (for the sake of parallel programming)
	std::vector<int> leafNodeList;
	for (const auto&[_nodeID, _parList] : _parMap){
		leafNodeList.push_back(_nodeID);
	}

	// Here we collect the data of target particle
	
	// Iterate through all node in the grid
	#pragma omp parallel for
	for (size_t n = 0; n < leafNodeList.size(); n++){
		// Here n is the node order in the leaf node list

		// Alias to the current node
		const int &_nodeID = leafNodeList[n];
		const Node* _currNode = _baseGrd.nodeMap.at(_nodeID);

		// Only evaluate the leaf node
        if(!(_currNode->isLeaf)) throw std::exception();
        
        // Get the list of neighbor node
        std::vector<Node*> _nghNodeList;
        _baseGrd.findNghAll(_nghNodeList, _currNode);

		// Aliasing to the current particle list
		const std::vector<int> &_parIDList = _parMap.at(_nodeID);

		// Iterate through all point inside the current node which the one for evaluation
		for (size_t i = 0; i < _parIDList.size(); i++){
			// Aliasing the current particle ID
			const int &ID_i = _parIDList[i];
			
            // Evaluate to each point inside the grid neighbor
			for (const auto &nghNode : _nghNodeList){
				for (const auto &ID_j : nghNode->parList){		
                    // Calculate the distance square between two points
					double _dr2 = 0;
					// Calculate in x direction
						double _dx = _par.x[ID_i] - _srcPar.x[ID_j];
						_dr2 += (_dx*_dx);
					// Calculate in y direction
						double _dy = _par.y[ID_i] - _srcPar.y[ID_j];
						_dr2 += (_dy*_dy);
					// Calculate in z direction
					if (DIM > 2){
						double _dz = _par.z[ID_i] - _srcPar.z[ID_j];
						_dr2 += (_dz*_dz);
					}

					// Evaluate distance
					if (Pars::opt_ngh_interact == 1){
						// Check the particle i
						double _rSup = _par.s[ID_i] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
						if (_dr2 < (_rSup*_rSup)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
					else if (Pars::opt_ngh_interact == 2){
						// Support radius of average size
						double _rSupAve = ((_par.s[ID_i] + _srcPar.s[ID_j]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
						if (_dr2 < (_rSupAve*_rSupAve)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
				}
			}
		}
	}
	return;
}

/**
 *  @brief  Evaluate the neighboring interaction for each particle. 
 *  Interaction pairs are determined by grid node relation.
 *  Given a single grid node that basically connect the target particle and
 *  take the relation of soruce particle toward the grid node using a map.
 *
 *  @param  _nghIDList	[OUTPUT] Neighbor ID list of eval particle toward source particle.
 *  @param  _trgPar 	The particle to evaluate neighbor. (The Target Particle for neighbor evaluation) [Eulerian Particle]
 *  @param  _srcPar		The source particle to evaluate neighbor. (Particle that hold the Properties value) [Lagrangian Particle]
 *  @param  _trgGridNode  	The node container of target particle. (The grid node that corresponding to target particle)
 *  @param  _srcParNodeMap  The source particle mapping to trgGridNode. (The node mapping of source particle to trgGridNode)
*/
void GridNodeNgh::eval_inter_ngh_gridNode(std::vector<std::vector<int>> &_nghIDList,
                                 const Particle &_trgPar,
                                 const Particle &_srcPar,
                                 const GridNode &_trgGridNode,
                                 const std::unordered_map<int, std::vector<int>> &_srcParNodeMap)
{
	// Justification of the particle size
	//  The size will directly using the target particle, 
	//  because the distribution between target and source is quite similar (source is only slightly advected)
	
	// **Evaluate neighbor interation
    _nghIDList.clear();
    _nghIDList.resize(_trgPar.num);

	// Collect all leaf node only (for the sake of parallel programming)
	std::vector<int> leafNodeList;
	for (auto &[_ID, _node] : _trgGridNode.nodeMap){
		if (_node->isLeaf){
			leafNodeList.push_back(_ID);
		}
	}
	
	// Iterate through all node in the grid
	#pragma omp parallel for
	for (size_t n = 0; n < leafNodeList.size(); n++){
		// Here n is the node order in the leaf node list

		// Alias to the current node
		const int &_nodeID = leafNodeList[n];
		const Node* _currNode = _trgGridNode.nodeMap.at(_nodeID);

		// Only evaluate the leaf node
        if(!(_currNode->isLeaf)) throw std::exception();
        
        // Get the list of neighbor node
        std::vector<Node*> _nghNodeList;
        _trgGridNode.findNghAll(_nghNodeList, _currNode);

		// Aliasing to the current particle list
		const std::vector<int> &_parIDList = _currNode->parList;

		// Iterate through all point inside the current node which the one for evaluation
		for (size_t i = 0; i < _parIDList.size(); i++){
			// Aliasing the current particle ID (Target Particle)
			const int &ID_i = _parIDList[i];
			
            // Evaluate to each point inside the grid neighbor
			for (const auto &nghNode : _nghNodeList){
				for (const auto &ID_j : _srcParNodeMap.at(nghNode->nodeID)){		
                    // Calculate the distance square between two points
					double _dr2 = 0;
					// Calculate in x direction
						double _dx = _trgPar.x[ID_i] - _srcPar.x[ID_j];
						_dr2 += (_dx*_dx);
					// Calculate in y direction
						double _dy = _trgPar.y[ID_i] - _srcPar.y[ID_j];
						_dr2 += (_dy*_dy);
					// Calculate in z direction
					#if DIM == 3
						double _dz = _trgPar.z[ID_i] - _srcPar.z[ID_j];
						_dr2 += (_dz*_dz);
					#endif

					// Evaluate distance
					if (Pars::opt_ngh_interact == 1){
						// Check the particle i
						double _rSup = _trgPar.s[ID_i] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
						if (_dr2 < (_rSup*_rSup)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
					else if (Pars::opt_ngh_interact == 2){
						// Support radius of average size
						double _rSupAve = ((_trgPar.s[ID_i] + _srcPar.s[ID_j]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
						if (_dr2 < (_rSupAve*_rSupAve)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
				}
			}
		}
	}
	return;
}


/**
 *  @brief  Assign the particle into the corresponding leaf node. This process
 *  evalaute the location of particle toward the grid node map, then take each node that
 *  correlated to the current evaluated particle. The map will groups all evaluated
 *  particles on the corresponding node.
 *
 *  @param  baseGrid The base node map data for node ID evaluation.
 *  @param  parNodeMap [OUTPUT] The map data of particle contained on node that consisted of leaf node only.
 *  This container maps the nodeID to all ID of evalulated particle inside the node 
 *  (or can be illustrated as : [nodeID]->{ [parID] }).
 *  @param  evalPar [UPDATE] Particle data to set the node ID.
*/
void GridNodeNgh::assign_par2node(const GridNode &nGrd, std::unordered_map<int, std::vector<int>> &parNodeMap, Particle &par){
    /** Procedure:
     *  1. Evaluate the ROOT node ID containing each particle
     *  2. Push down the node ID until it reaches the leaf node
    */

    // Node container variable
    std::vector<int> nodeList1, nodeList2;          // The temporary node ID container
    std::vector<int> *currNodeList, *nextNodeList;  // The temporary alias of node ID container

	// Reserve the particle node map container first
	parNodeMap.clear();

    // PROCEDURE 1:
    // ***********
    // Determine the ID of the corresponding root node
    par.nodeID.resize(par.num, 0);  // Reserve the node ID container
    for (int parID = 0; parID < par.num; parID++){
        // Retrieve the particle coordinate
        double coord[DIM];
        coord[0] = par.x[parID];
        coord[1] = par.y[parID];
        #if (DIM == 3 )
        coord[2] = par.z[parID];
        #endif

        // Get the root ID
        int _nodeID = nGrd.pos2ID(coord, ROOT_LEVEL);
        par.nodeID[parID] = _nodeID;

        // Update the temporary container
        if (parNodeMap.count(_nodeID) == 0){
            nodeList1.push_back(_nodeID);
        }
        
        // Update the nodeID to parID map container data
        parNodeMap[_nodeID].push_back(parID);
    }


    // PROCEDURE 2:
    // ***********
    // Iteratively evaluate leaf and downpass to child for non leaf node

    // Start iteration
    for (int level = 0; level < Pars::max_level; level++){
        // Create aliasing for the node list
        if (level % 2 == 0){
            // Even level
            currNodeList = &nodeList1;
            nextNodeList = &nodeList2;
        }else{
            // Odd Level
            currNodeList = &nodeList2;
            nextNodeList = &nodeList1;
        }

        // Set up the necessary container
        nextNodeList->clear();

        // Check each node in the current list
        for (const int &nodeID : *currNodeList){            
            // Check the current node is not a leaf
            if (nGrd.nodeMap.at(nodeID)->isLeaf != true){
                // Distribute the particle into the correspond child (if current node is NOT a leaf)

                // Internal container
                std::vector<int> &parList = parNodeMap.at(nodeID);      // List of all particle in the current node
                std::vector<int> chdIDList = nGrd.findChild(nodeID);    // List of all child node ID of the current node

                // Coordinate container
                double nodeCoord[DIM];  // Node mid point coordinate
                double coord[DIM];      // Particle coordinate

                // Retrieve the middle coordinate of current node
                const Node &_node = *(nGrd.nodeMap.at(nodeID));
                basis_loop(d) nodeCoord[d] = _node.pivCoor[d] + (0.5 * _node.length);
                
                // Evaluate all particle in the current node
                for (const int &parID : parList){
                    // Retrieve the particle coordinate
                    coord[0] = par.x[parID];
                    coord[1] = par.y[parID];
                    #if (DIM == 3 )
                    coord[2] = par.z[parID];
                    #endif

                    /** Illustration of LID (local ID) numbering
                     *   2D space       y                                    *  Comparing value in the binary means!
                     *   _________     |                                     *  -----------------------------------
                     *  | 2  | 3  |    |____ x                               *    2D space only have 2 basis x,y
                     *  |____|____|                    Back Block            *    it can be translated to 2 digit binary
                     *  | 0  | 1  |    3D space      ___________________     *     For clarity, apply the LID number in 2D space
                     *  |____|____|                /:        /:        /|    *     0 -> 00 : x left to center  | y below center
                     *                            /_:______ /_:______ / |    *     1 -> 01 : x right to center | y below center
                     *       Front Block         |  :  2   |  :  3   |  |    *     2 -> 10 : x left to center  | y above center
                     *      ___________________  |  :......|..:......|..|    *     3 -> 11 : x right to center | y above center
                     *    /:        /:        /| | ,:      | ,:      | /|    *          ^^
                     *   /_:______ /_:______ / | |,_:______|,_:______|/ |    *       y--||--x
                     *  |  :  6   |  :  7   |  | |  :  0   |  :  1   |  |    *  
                     *  |  :......|..:......|..| |  :......|..:......|..|    *  The first digit (rightmost digit) represents
                     *  | ,:      | ,:      | /| | ,       | ,       | /     *   the relative x location of child toward parent node center, 
                     *  |,_:______|,_:______|/ | |,________|,________|/      *   otherwise the second digit represent the relative y location.
                     *  |  :  4   |  :  5   |  |        y                    *  
                     *  |  :......|..:......|..|       |                     *  The translation value between child relative position
                     *  | ,       | ,       | /        |____ x               *   and the LID value is used to determine the child node ID
                     *  |,________|,________|/        /                      *  
                     *                             z /                       *  This idea works fine for the extension to 3D!  
                    */
                    int LID = 0;    // The local ID position of child location [initialized by 0]

                    // Determine the child local position [quad (2D) or octant (3D)] by
                    // comparing the particle coordinate toward current node mid point
                    basis_loop(d){
                        LID |= ((coord[d] > nodeCoord[d] ? 1 : 0) << d);
                        // NOTE:
                        // > Use the binary left shift (<<) to set the digit location
                        // > Use the binary or operator to overwrite binary value to LID
                    }

                    // Get the new ID
                    int newNodeID = chdIDList[LID];     // The new node ID
                    par.nodeID[parID] = newNodeID;      // Update the node ID of current particle

                    // Update the temporary container
                    if (parNodeMap.count(newNodeID) == 0){
                        nextNodeList->push_back(newNodeID);
                    }

                    // Update the nodeID to parID map container data
                    parNodeMap[newNodeID].push_back(parID);
                }

                // Release the container at current node
                parNodeMap.erase(nodeID);
            }
        }
    }

    /** Procedure:
     *  1. Put all particle into the corresponding root node (also take the other one)
     *     -> For leaf node     : Keep as settled particle
     *     -> For non leaf node : Transfer to the corresponding child
     *  2. Do the first and second check
     *  3. Done when the level reach max (automatically a leaf node)
    */

    return;
}