#include "gridNodeAdapt.hpp"
// #include "../save_data/save_data.hpp"



// ===================================================== //
// +---------------- Adaptation Method ----------------+ //
// ===================================================== //
// #pragma region ADAPTATION_METHOD

/**
 *  @brief  Perform the node and particle adaptation based on the value tolerance evaluation.
 *         
 *  @param  baseGrid [UPDATE] The grid to perform adaptation operation.
 *  @param  tempPar  The eulerian particle to hold the basic particle position and particle relation (neighborhood). 
 *  @param  currPar  The lagrangian particle that contains the physical property.
 *  @param  nghConv  Perform neighbor data convertion to the calculation. [0:= Evaluate neighbor, other not]
 *  @param  errPredProp  The property value for adaptation evaluation.
 *  @param  parNodeMap   The mapping of currPar toward the baseGrid node.
 *  
 *  @return The adaptation flag, said that adaptation is performed.
*/
void GridNodeAdapt::get_adaptation_LSMPS(GridNode &baseGrid, const Particle &currPar, 
                              std::vector<std::vector<double>*>& errPredProp){
    // Initial Adaptation Log
    MESSAGE_LOG << "Initial total nodes    : " << baseGrid.nodeMap.size() << "\n";
    MESSAGE_LOG << "Initial total particle : " << currPar.num << "\n";

    /* PROCEDURE !!
        1. Set up and initialize the parameter in adaptation procedure
        2. Refine all nodes in refinement list
        3. Evaluate NLD of all nodes near refined node
        4. Compress all nodes in compress list
        5. Update the particle and GridNode data [Update the other data]

        NOTE:
        > At PROCEDURE [2], [3], & [4] the new particle data and 
          neighbor (toward old particle) is generated
        > The temporary grid will be deleted in the very end
    */

    // PROCEDURE 1!
    // ************
    // Set up and initialize the parameter in adaptation procedure
    
    // Time counter
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // // Set the selected property to be evaluated for adaptation
    // // *NOTE : Change the property if necessary
    // std::vector<double> selProp = currPar.vorticity;
    // for (size_t i = 0; i < selProp.size(); i++){
    //     selProp[i] = std::abs(selProp[i]);
    // }

    // // EVALUATION 1
    // // Find the maximum value of the selected properties
    // this->maxValue = *(std::max_element(selProp.begin(),selProp.end()));
    // // MESSAGE_LOG << "The selected property max. value : " << this->maxValue << "\n\n";

    // EVALUATION 2
    // MESSAGE_LOG << "Set target level of each node in grid!\n";
    // Set the target resolution level of each node in baseGrid
    this->targetLevel(baseGrid, currPar, errPredProp);

    // EVALUATION 3
    // MESSAGE_LOG << "Set adaptation flag!\n";
    // Set adaptation flag and group each LEAF leaf node
    this->set_adaptation_flag(baseGrid);

    // Time counter
    // Display calculation time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<+> Adapt 1: Initialize    : %f s\n", _time);

    // // [DEBUG] Save the temporary grid
    // baseGrid.saveLeafGrid(tempGrid, "Adaptation");

    // // ========================================
    // // ------------ FIRST DEBUGGER ------------
    // // ========================================
    // std::cout << "\n[LOG 1] PRINT NUMBER OF NODE!\n";
    // std::cout << FONT_GRAY;
    // std::cout << " >>> Refinement  : " << this->refineList.size() << "\n";
    // std::cout << " >>> IDLE        : " << this->idleList.size() << "\n";
    // std::cout << " >>> Compression : " << this->compressList.size() << "\n\n";
    // std::cout << FONT_RESET;

    // No need adaptation if there is no refinement and compression
    if (this->refineList.empty() && this->compressList.empty()){
        WARNING_LOG << "Be carefull, adaptation is cancelled!!!\n";
        return;
    }

    // ================ ======================
    // +----- PARTICLE INITIALIZATION ------+
    // ======================================
    // [NOTE] : Modify this line if necessary for next code refactor
    // Create a new particle and the alias
    // this->newPar = new Particle;
    // Particle &_newPar = *this->newPar;

    // // Initialize all data in new particle
    // _newPar.num = 0;
    // _newPar.x.clear();
    // _newPar.y.clear();
    // _newPar.z.clear();
    // _newPar.nodeID.clear();
    // _newPar.neighbor.clear();
    // _newPar.level.clear();
    // _newPar.s.clear();         // The only variable that follow the predecesor
    // _newPar.isActive.clear();  // The second variable that follow the predecesor

    // CLASS GLOBAL PARAMETER
    // **********************
    // This must be put into the class member
    // Initialization data for particle generation (USED IN REFINEMENT, NLD, AND COMPRESSION)
    this->parNum = baseGrid.baseParNum;         // The number of particle in a node in one dimension
    this->totParNum = Pars::intPow(this->parNum, DIM);      // The total number of particle in a node

    // Initialize the divisor
    this->div[0] = 1;
    for (int i = 1; i < DIM; i++){
        this->div[i] = this->div[i-1] * this->parNum;
    }

    /** NARRATION:
     *   In this adaptation proces
     * 
     *  Before Adaptation (Grid Node)             After Adaptation (Grid Node)             
     *   ___________ ___________ ___________       ___________ ___________ ___________      
     *  |           |+          |-    |-    |     |           |*    |*    |#          |     
     *  |           |           |_____|_____|     |           |_____|_____|           |     
     *  |           |           |-    |-    |     |           |*    |*    |           |     
     *  |___________|___________|_____|_____|     |___________|_____|_____|___________|     
     *  |     |     |+    |+    |     |     |     |     |     |*_|*_|*_|*_|     |     |     
     *  |_____|_____|_____|_____|_____|_____|     |_____|_____|*_|*_|*_|*_|_____|_____|     
     *  |     |     |     |     |     |     |     |     |     |     |     |     |     |     
     *  |_____|_____|_____|_____|_____|_____|     |_____|_____|_____|_____|_____|_____|     
     * 
     *  Description:
     *   > Before the adaptation Node with * is refined node, while the one with # is compressed node
     *      -> Node with + will be refined, 
     *      -> A group of sibling node with - will be compressed
     *   > After the adaptation
     *      -> Node with * is refined node, 
     *      -> Node with # is compressed node
     * 
     *  In order to process the adaptation
     *   > Refinement will be processes in different node called @param TemporaryGrid
     *   > The Compression will be done directly on the @param BaseGrid
     *   > The particle data will be generated accordingly in @param newParticle
     * 
     *  @param Refinement -> All performed in the @param TemporaryGrid
     *   +) The process will be done recursively for each node, e.g.
     *           ___________          ___________          ___________ 
     *          |+          |        |     |     |        |     |     |
     *          |           |   ->   |_____|_____|   ->   |_____|_____|
     *          |           |        |     |+    |        |     |__|__|
     *          |___________|        |_____|_____|        |_____|__|__|
     *            Head Node             Refine 1             Refine 2
     *        *) Only refine the node montioned by (+) sign
     * 
     *   +) Take head node -> find all leaf neighbor and put into @param IDqueue for NLD evaluation
     *   +) Refine head node until meet target level (RECURSIVELY; each refinement need to check TRL)
     *   +) Generate the @param newParticle at each leaf of refined node
     * 
     *  @param NLD (Neighbor Level Different Evaluation)
     *   +) Illustration of level window (e.g. Max Level = 4, NLD = 1)                                              
     *          @param ROOT_LEVEL    0     ^
     *            ------------             |  In this level bound, we need to 
     *                 :             1     |   evaluate the NLD criteria
     *                 :             2     v
     *            ------------            
     *          @param NLD_BOUND     3    --- From here on we don't need to
     *            ------------             |   evaluate the NLD criteria
     *          @param MAX_LEVEL     4     v  
     *            ------------
     *   +) Must evaluate from (Max level - 2)
     *   +) At each node, find neighbor node, 
     *       -> neighbors with lower level put into next NLD evaluation,
     *       -> neighbors with higher level as current NLD eval (refine or not)
     *   +) If need refine (see the refinement process -> child will be put back into NLD evaluation)
     *   *) PS: @param ParticleDivision method is used to move the source particle into refined child node
     * 
     *  @param Compression
     *   +) Sort the node, only collection of sibling node that need compression will be compressed
     *   +) Check the parent NLD, if meet the criteria, the group of node is allowed to be compressed
     *   +) Compress the node and generate particle
     * 
    */


    // PROCEDURE 2!
    // ************
    // Refine all nodes in refinement list, for each node (called a HEAD node):
    // [1] Check the LEAF node neighbor of the HEAD node
    // [2] Refine the Head node
    // [3] Evaluate refinement for all child until the LEAF node (A queue evaluation)
    //      -> Find the TRL of current node
    //      -> Refine the node if meet criteria
    //      -> Generate particle if meet criteria

    // ID container of [QUEUE LIST] (*For NLD criteria evaluation in PROCEDURE 3)
    const int LEVEL_BOUND = this->maxLevel - this->NghlevelDiff;    // The upper bound of level limit to be evaluated in queue (*not including the bound)
    std::vector<std::vector<int>> IDqueue(LEVEL_BOUND);             // For NLD evaluation
    std::unordered_map<int,bool> IDflag;    // Flag is rather [true] := in base grid; or [false] := in temporary grid

    // Refinement ID container
    std::vector<int> nodeIDrefine;          // A queue list of all node ID that are to be refined
    
    // START PROCEDURE 2:
    int refinedNode = 0;        // Refinement counter for console log

    // Time counter
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif
    
    // MESSAGE_LOG << "Grid node refinement by adaptation!\n";
    int iter = 0;               // For debugging saving file
    
    // [LOOP 1] -> Loop through all node in refine List
    for (auto &[_headID,_headFlag] : this->refineList){
        // Internal container
        std::vector<int> tempChdIDList;     // Temporary child ID list container

        // Create an alias to the current node as a HEAD NODE
        Node *&_headNode = this->tempGrid.nodeMap.at(_headID);  // This node is indeed a LEAF

        // // An update toward the particle division
        // _headNode->parList = parNodeMap.at(_headNode->nodeID);
        
        // **[P.1] Check each neighbor
        // *******
        //    -> [A] Evaluate all node neighbor (on the baseGrid) of the HEAD node
        //    -> [B] At the same time, collect all particle inside all neighbor node
        // Create a particle neighbor candidate ID list
        std::vector<Node*> nodeNghIDList;   // The neighbor list (consisted of LEAF node) at "baseGrid"
        // std::vector<int> parNghIDList;      // List of all particle inside the neighbor (collect the lagrangian particle)
        
        // * Find all LEAF Node neighbor
        baseGrid.findNghAll(nodeNghIDList, _headNode);

        // * Check all neighbor node (For NLD evaluation)
        for (auto &_nghNode : nodeNghIDList){
            // *[A] Check the neighbor if it meets the NLD criteria
            int &_nodeNghID = _nghNode->nodeID;
            // Check the neighbor level
            if (_nghNode->level < LEVEL_BOUND){
                // Only evaluate the node outside refinement list
                if (this->refineList.count(_nodeNghID) == 0){
                    // The node is outside REFINARY ID list container
                    // Check whether the node is already existed in the list queue
                    if (IDflag.count(_nodeNghID) == 0){
                        // Only put into the queue if the node is not existed
                        IDqueue[_nghNode->level].push_back(_nodeNghID);
                        IDflag.insert({_nodeNghID, true});      // TRUE flag, the node inside "baseGrid"
                    }
                }
            }

            // // *[B] Push all particle inside the current neighbor node into the particle ID list
            // // MUST BE UPDATED
            // for (auto &_nghParID : _nghNode->parList){
            //     parNghIDList.push_back(_nghParID);
            // }
            // if (parNodeMap.count(_nodeNghID) != 0){     // Make sure the lagrangian particle is inside the node
            //     for (auto &_nghParID : parNodeMap.at(_nodeNghID)){
            //         parNghIDList.push_back(_nghParID);
            //     }
            // }
        }

        // **[P.2] Refine the head node
        // *******
        tempChdIDList.clear();
        this->tempGrid.refineNode(_headNode, baseGrid, tempChdIDList);      // Refine the head node (based on the baseGrid)
        this->tempGrid.divideParticle(_headNode, currPar, tempChdIDList, 1);    // Divide the particle from head node to child node

        // !!!! Assign particle to evaluate the TRL !!!!

        refinedNode++;      // Node is refined add into account

        // * Put the first child into the queue list
        // Reserve the space of refine node ID container
        nodeIDrefine.clear();       // A recursive loops (until there is no node to be refined)
        // Put the child into refinement list (A queue list)
        for (auto &chdID : tempChdIDList) nodeIDrefine.push_back(chdID);

        // **[P.3] Evaluate refinement for all child (RECURSIVE LOOP)
        // *******
        // [LOOP 2] -> Loop through all child to be evaluated for refinement
        for (size_t i = 0; i < nodeIDrefine.size(); i++){
            // Procedure in evaluating each child
            //    [A] Set the target refinement level of the node (use the particle in parList)
            //    [B] Check the refinement flag
            //        [COND 1] REF the eval Node that TLR is higher than the current level
            //        [COND 2] Put into NLD queue if the eval node level is not at maximum level**
            //        [COND 3] PG on the node that is not going into the NLD check
            // *Desc. [REF]->Refine node; [PG]->Particle generation and neighbor evaluation
            // **The maximum level is the MAX_NGH_DIFF_LEVEL below than the maximum level 
            
            // Create alias to the current node ID
            int &_evalID = nodeIDrefine[i];
            // A hanlDer to prevent further problem
            if (this->tempGrid.nodeMap.count(_evalID) == 0){
                ERROR_LOG << "Node " << _evalID << " is not existed!\n";
                throw std::runtime_error("ERROR [ADAPTATION]#1 : The evaluated node is not existed!");
            }
            // Create alias to the current node
            Node *&_evalNode = this->tempGrid.nodeMap.at(_evalID);

            // *[A] Check the target level TRL of the current evaluated Node 
            // Check from all particles in the node particle list
            for (int &_particleID : _evalNode->parList){
                // int tarLv;
                // double const &currVal = selProp[_particleID];
                // double tolerance = this->maxValue * Pars::adapt_tol;
                // for (tarLv = baseGrid.maxLevel; tarLv > 0; tarLv--){
                //     if (currVal > tolerance){
                //         break;
                //     }else{
                //         tolerance *= Pars::adapt_tol;
                //     }
                // }

                // // Target level of the current evaluated node is the maximum TRL between all particles inside the node
                // _evalNode->tarResLv = std::max<int>(_evalNode->tarResLv, tarLv);

                this->targetLevelCalc(_evalNode, _particleID, errPredProp);
            }
            
            
            // *[B] Check the refinement flag of the current node
            // [COND 1] The current level is lower than target level
            if(_evalNode->level < _evalNode->tarResLv){
                // Refine the current evaluated node
                tempChdIDList.clear();
                this->tempGrid.refineNode(_evalNode, baseGrid, tempChdIDList);
                this->tempGrid.divideParticle(_evalNode, currPar, tempChdIDList, 1);

                // Need an assign particle for TRL
                // Divide Particle is utilized for refinement not for neighbor evaluation (But then check the compression)

                refinedNode++;      // Node is refined add into account
                
                // Put the child into refinement list (A queue list)
                for (auto &__chdID : tempChdIDList) nodeIDrefine.push_back(__chdID);

            }
            // [COND 2] The level of current node still in range of NLD evaluation
            else if (_evalNode->level < LEVEL_BOUND){
                // Error hanlDer
                if (IDflag.count(_evalNode->nodeID) != 0){
                    ERROR_LOG << "Node " << _evalNode->nodeID << " is already in the list!\n";
                    throw std::runtime_error("ERROR [ADAPTATION]#2 : The current node is already in the queue list!");
                }
                
                // Put the current node ID into the NLD queue list
                IDqueue[_evalNode->level].push_back(_evalNode->nodeID);
                IDflag.insert({_evalNode->nodeID, false});    // FALSE flag, the node inside "temporary Grid"
            }
            // [COND 3] Particle generation for the node outside the previous criteria
            else{
                // // Generate the particle inside the current node
                // double _headSize = _headNode->length / this->parNum;
                // this->generateParticle(_evalNode, _headSize);

                // // Evaluate the neighbor
                // if (nghConv == 0)
                // this->findParticleNeighbor(_evalNode->parList, parNghIDList, currPar);
            }
        }
    }

    // Display calculation time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<+> Adapt 2: Refinement    : %f s\n", _time);
    
    
    // // [DEBUG] Check the temporary grid
    // baseGrid.saveLeafGrid(tempGrid,"AdaptationRefined");
    // save_particle(*this->newPar, "AdaptationRefined");

    // // ========================================
    // // ----------- SECOND DEBUGGER ------------
    // // ========================================
    // std::cout << "\n[LOG 2] PRINT NUMBER OF NODE!\n";
    // std::cout << FONT_GRAY;
    // std::cout << " >>> Refinement  : " << this->refineList.size() << "\n";
    // std::cout << " >>> IDLE        : " << this->idleList.size() << "\n";
    // std::cout << " >>> Compression : " << this->compressList.size() << "\n\n";
    // std::cout << FONT_RESET;


    // PROCEDURE 3!
    // ************
    // Neighbor Level Different criteria evaluation 
    //   *see the PROCEDURE 4 at "generateGrid.hpp"
    
    // Refine the node according to neighbor level different (NLD) criterion
    // Do a loop leveled queue check through all nodes that need an NLD evaluation
    
    // Procedure Operation:
    // [1] Loop check until there are no nodes left on the queue
    //        > Queue container [IDqueue] and [IDflag]
    // [2] Loop through all nodes on each level (from highest level in the queue)
    // [3] At each node check the neighbor node ID (Further operation see inside the code)

    // Debug FLAG
    bool DEBUG = false;
    iter = 0;
    
    // START PROCEDURE 3:
    // MESSAGE_LOG << "Evaluate NLD criteria and refinement\n";
    
    // Time counter
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif
    
    // LOOP [1]: Loop until no node left on the queue
    bool stop = false;      // Loop trigger
    int fineLvl = this->maxLevel - this->NghlevelDiff - 1;  // Starting level at each loop
    while(!stop){
        if (DEBUG) std::cout << FONT_GREEN << "[LOG] Outer while Loop " << iter << "\n" << FONT_RESET;
        
        // Activate the trigger of loop termination
        // The loop need to halt when the queue is empty
        stop = true;
        
        // LOOP [2]: Check the ID inside the queue from the finest level
        for (int currLvl = fineLvl; currLvl >= 0; currLvl--){
            if (DEBUG) std::cout << FONT_GREEN << "[LOG] Iteration for the level " << currLvl << "!\n" << FONT_RESET;
            
            // Cancel the termination trigger if the queue is not empty
            if (!(IDqueue[currLvl].empty())){
                stop = false;
            }

            // The "IDqueue" container consisted of 2 flag:
            // [TRUE]   -> leaf node from compresed and idle list (baseGrid)
            // [FALSE]  -> leaf node from refine ID list (tempGrid)

            // LOOP [3]: Check the Node inside the queue at current iteration level
            // We can do parallel here
            for (int &_ID : IDqueue[currLvl]){
                // Node Evaluation Procedure
                // [0] Check if leaf node : only evaluate leaf node
                // [1] Find the neighbor ID at the current level
                // [2] Check the node at neighbor ID
                //        > [COND 1] The finest neighbor level different is larger than MAX_LEVEL_DIFF
                //            -> Refine the current node -> Put the refined node into container
                //        > [COND 2] Existed: The neigbor with level lower than current level
                //            -> Put the ID into the queue in accordance to its level

                // Create a node alias
                Node *_Node;
                // Find the location of _Node
                if(IDflag.at(_ID)){ 
                    _Node = baseGrid.nodeMap.at(_ID);
                    // *Need to refine or generate particle also transfer the list
                }else{
                    _Node = this->tempGrid.nodeMap.at(_ID);
                    // *Only need to refine or generate particle
                }

                // ** [0] Check leaf node <?> Not needed, already fulfilled <?>
                if (!(_Node->isLeaf)){
                    // Proceed to the next Node in queue
                    WARNING_LOG << " NLD check on a non-LEAF node! [CODE:0xNG!!]\n";
                    continue;
                }

                // ** [1] Find the neighbor ID at the current level
                std::vector<int> nghIDList;             // List of the neighbor ID
                std::unordered_map<int,int> nghIDflag;  // Flag of existed neighbor [ID, flag]
                                                        //   -1:= node located at tempGrid
                                                        //    0:= node is not existed
                                                        //    1:= node located at baseGrid
                
                // Update each list
                baseGrid.findNghLvl(nghIDList, _Node);  // Update neigbor list
                for (int &_nghID : nghIDList){          // Update neigbor flag
                    if (this->tempGrid.nodeMap.count(_nghID) != 0){
                        nghIDflag.insert({_nghID, -1});
                        // A node that have been refined [type 1] or the node head at refined [type 2]
                    }else if (baseGrid.nodeMap.count(_nghID) != 0){
                        nghIDflag.insert({_nghID, 1});
                    }else{
                        nghIDflag.insert({_nghID, 0});
                    }
                }

                // ** [2] Check each neigbor node (Find the leaf node + check need refinement or not)
                bool needRefine = false;
                // LOOP [4]: Check each neigbor node
                for (size_t i = 0; i < nghIDList.size(); i++){
                    // Retrieve the neighbor ID
                    int _nghID = nghIDList[i];

                    // No need to evaluate to itself
                    if (_nghID == _ID) continue;

                    // Evaluation exeption!!!
                    if (needRefine){
                        if(_nghID >= baseGrid.getPivID(_Node->level+1)){
                            // Skip the evaluation for finer neighbor
                            //    if the current node is need to be refined
                            // If the current node flag to refine have been settled
                            //    it no need to check higher level Node
                            continue;
                            // Continue so we can evaluate the lower level neighbor node to put into next NLD evaluation
                        }
                    }
                    
                    // Evaluation:
                    // [A] Node not existed
                    //       -> Navigate to the parent node
                    // [B] Node existed
                    //       [a] If leaf node
                    //           -> Evaluate the node [COND 1] & [COND 2]
                    //       [b] If not a leaf node
                    //           -> Navigate to child

                    // ** [A] Check whether the node is not existed
                    if (nghIDflag.at(_nghID) == 0){ // Node is NOT existing
                        // To cancel recursive calculation
                        // Check whether the candidate node level is bigger than the current level
                        if (_nghID >= baseGrid.getPivID(currLvl + 1)){
                            // There is no way a Node is not existed when its ID went to a higher level
                            WARNING_LOG << " An attempt UP to a recursive check is almost happenned!\n";
                            continue;
                        }
                        
                        // Find parent ID
                        int nghParID = baseGrid.findParent(_nghID);
                        
                        // Put the parent node ID into the neighbor ID list to be evaluated further
                        if (nghIDflag.count(nghParID) == 0){
                            // Only put the node when it is not existing in the list
                            nghIDList.push_back(nghParID);
                            
                            // Check whether it is located at temporary or base grid
                            if (baseGrid.nodeMap.count(nghParID) != 0){
                                // Parent node of neighbor lies on baseGrid
                                nghIDflag.insert({nghParID, 1});
                            }else if (this->tempGrid.nodeMap.count(nghParID) != 0){
                                // Parent node of neighbor lies on tempGrid
                                nghIDflag.insert({nghParID, -1});
                                // // [!] It is impossible to be located in this region
                                // ERROR_LOG << "An Impossible data structure is occured on node " << _Node->nodeID << "\n";
                                // throw std::runtime_error("ERROR [NDL EVAL] Targeting on impossible structure!");
                                // Anyway this would be able when the surrounding is refined node
                            }else{
                                // Parent node of neighbor is not existed
                                nghIDflag.insert({nghParID, 0});
                            }
                        }
                    }
                    
                    // ** [B] Node is existing, Check whether a leaf node (*at base grid)
                    else if (nghIDflag.at(_nghID) == 1){ // This line must be safe, because it will check the first condition first, only check second if first is true
                        if (baseGrid.nodeMap.at(_nghID)->isLeaf){ // Node is existing in the baseGrid
                            // Check the level of the current neighbor node
                            int _nghLvl = baseGrid.getLevel(_nghID);

                            // ** [COND 2] Put the ID into queue container (*if meet the criteria)
                            // Whatever the node is, if the level is lower we need to check on next NLD
                            if (_nghLvl < _Node->level){
                                // Add the node ID to the queuqe
                                if (IDflag.count(_nghID) == 0){
                                    // Only add if the ID is not existing in the queue
                                    IDqueue[_nghLvl].push_back(_nghID);
                                    IDflag.insert({_nghID, true});
                                }
                            }
                        }else{
                            // This section is imposible
                        }
                    }
                    
                    // ** [B] Node is existing, Check whether a leaf node (*at temporary grid)
                    else if (nghIDflag.at(_nghID) == -1){ // This line must be safe, because it will check the first condition first, only check second if first is true
                        if(this->tempGrid.nodeMap.at(_nghID)->isLeaf){ // Node is existing in the temporary grid
                            // Check the level of the current neighbor node
                            int _nghLvl = baseGrid.getLevel(_nghID);

                            // ** [COND 2] Put the ID into queue container (*if meet the criteria)
                            if (_nghLvl < _Node->level){
                                // Add the node ID to the queuqe
                                if (IDflag.count(_nghID) == 0){
                                    // Only add if the ID is not existing in the queue
                                    IDqueue[_nghLvl].push_back(_nghID);
                                    IDflag.insert({_nghID, false});
                                }
                            }
                        }else{
                            // Check the child
                        }
                             
                    }
                    
                    // else if (baseGrid.nodeMap.at(_nghID)->isLeaf){ // Node is existing in the baseGrid
                    //     // Check the level of the current neighbor node
                    //     int _nghLvl = baseGrid.getLevel(_nghID);

                    //     // ** [COND 2] Put the ID into queue container (*if meet the criteria)
                    //     // Whatever the node is, if the level is lower we need to check on next NLD
                    //     if (_nghLvl < _Node->level){
                    //         // Add the node ID to the queuqe
                    //         if (IDflag.count(_nghID) == 0){
                    //             // Only add if the ID is not existing in the queue
                    //             IDqueue[_nghLvl].push_back(_nghID);
                    //             IDflag.insert({_nghID, true});
                    //         }
                    //     }
                    // }
                    
                    // ** [C] Proceed to check child if the current node still no need to refine
                    else if (!needRefine){
                        // To cancel recursive calculation
                        // Check whether the candidate node level is smaller than the current level
                        if (_nghID < baseGrid.getPivID(currLvl)){
                            // There is no way a Node is not a leaf when existed and have no children
                            WARNING_LOG << " An attempt DOWN to a recursive check is almost happenned!\n";
                            continue;
                        }
                        
                        // Evaluate whether the current node need to be refined
                        int chdLvl = baseGrid.getLevel(_nghID) + 1;
                        if (chdLvl > _Node->level + this->NghlevelDiff){
                            // Meet the criteria to refine
                            needRefine = true;  // Here is the refine for current node, not the negihbor node
                            continue;           // Continue evaluate next neighbor node
                        }

                        // Find all children ID
                        std::vector<int> chdIDList;
                        if(nghIDflag.at(_nghID) == 1){
                            chdIDList = baseGrid.findChild(baseGrid.nodeMap.at(_nghID));
                        }
                        else if(nghIDflag.at(_nghID) == -1){
                            chdIDList = baseGrid.findChild(this->tempGrid.nodeMap.at(_nghID));
                        }
                        
                        // Put the child IDs into the candidate neighbor queue container (to be evaluated further)
                        for(int &chdID : chdIDList){
                            // Put the node to the candidate neighbor container
                            nghIDList.push_back(chdID);

                            // Check where is the head located at (temporary or base Grid)
                            if(nghIDflag.at(_nghID) == 1){
                                // Check whether the child inside refine list
                                if(this->refineList.count(chdID) == 0){ // Not at the refined List
                                    // Child located at baseGrid
                                    nghIDflag.insert({chdID, 1});
                                }else{
                                    // Child located at tempGrid
                                    nghIDflag.insert({chdID, -1}); // Flag changed from 1 to -1
                                }
                            }
                            else if(nghIDflag.at(_nghID) == -1){
                                // Child located at tempGrid
                                nghIDflag.insert({chdID, -1});
                            }
                        }
                    }

                } /* LOOP [4]: Check each neigbor node */

                
                // Refinement procedure in NLD:
                //    [1] TRUE flag   : Located at compression/idle list
                //                        Need to change the 
                //    [2] FALSE flag  : Located at refine list

                // ** [COND 1] Perform the refinement (*if meet the criteria)
                if (needRefine){
                    // Additional modification for node in baseGrid (need to move in refinement list and temporary grid)
                    if(IDflag.at(_ID)){
                        // Insert the current node into the refine ID list
                        this->refineList.insert({_ID,true});

                        // * Remove the current node from the compression / idle list
                        // At the compression list
                        if (this->compressList.count(_ID) != 0){
                            baseGrid.nodeMap.at(_ID)->needCompression = false;  // Turn off the flag
                            this->compressList.erase(_ID);
                        }
                        // At the idle list
                        else if (this->idleList.count(_ID) != 0) this->idleList.erase(_ID);

                        // Create a duplicate node to the current node
                        Node *_dupNode = new Node(_Node);
                        _Node = _dupNode;   // Change the node reference to the duplicate node
                        _Node->headNodeID = _Node->nodeID;  // Put the current ID as the head node
                        
                        // Push the duplicate node into the temporary grid
                        this->tempGrid.nodeMap.insert({_ID,_Node});
                    }
                    
                    // Refine the current node
                    std::vector<int> chdIDList;     // Container to store the child ID
                    this->tempGrid.refineNode(_Node, baseGrid, chdIDList);
                    // Update the head node of current node ID and all child
                    for (int &chdID : chdIDList){
                        Node*& chdNode = this->tempGrid.nodeMap.at(chdID);
                        chdNode->headNodeID = _Node->headNodeID;
                    }
                    _Node->headNodeID = -1;
                    // this->tempGrid.divideParticle(_Node, tempPar, chdIDList, 1);
                    refinedNode++;      // Node is refined add into account

                    // Store the child ID into the queue container (*if inside the NLD evaluation criteria)
                    if (currLvl + 1 < LEVEL_BOUND){
                        // Must be not existed inside the queue yet because the node was just created
                        for (int &chdID : chdIDList){
                            IDqueue[currLvl+1].push_back(chdID);
                            IDflag.insert({chdID, false});  // Inside the temporary list
                        }
                    }
                    
                    // Generate the particle inside each child node (*if not meet the criteria)
                    else{
                        // Find all particle at the neighbor nodes of the "HEAD node"

                        // Aliasing the node
                        Node *&_pivChdNode = this->tempGrid.nodeMap.at(chdIDList[0]);       // Aliasing the pivot child ID node
                        Node *&_headNode = baseGrid.nodeMap.at(_pivChdNode->headNodeID);    // Aliasing the HEAD node
                        
                        // // Create a particle neighbor candidate ID list
                        // std::vector<Node*> nodeNghIDList;
                        // std::vector<int> parNghIDList;
                        
                        // // Find all LEAF neighbor at the HEAD node
                        // baseGrid.findNghAll(nodeNghIDList, _headNode);
                        // for (auto &_nghNode : nodeNghIDList){
                        //     // // Push all particle into the neighbor particle ID list
                        //     // for (auto &_nghParticleID : _nghNode->parList){
                        //     //     parNghIDList.push_back(_nghParticleID);
                        //     // }
                        //     if (parNodeMap.count(_nghNode->nodeID) != 0){     // Make sure the lagrangian particle is inside the node
                        //         for (auto &_nghParID : parNodeMap.at(_nghNode->nodeID)){
                        //             parNghIDList.push_back(_nghParID);
                        //         }
                        //     }
                        // }
                        
                        // // Iterate through all child Node (generate the particle and neighbor)
                        // for (int &_chdID : chdIDList){
                        //     // Aliasing the current node
                        //     Node *&_currNode = this->tempGrid.nodeMap.at(_chdID);

                        //     // Generate the particle
                        //     double _headSize = _headNode->length / this->parNum;
                        //     this->generateParticle(_currNode, _headSize);
                            
                        //     // Find the neighbor
                        //     if (nghConv == 0)
                        //     this->findParticleNeighbor(_currNode->parList, parNghIDList, currPar);
                        // }
                    }
                }

                // // ** Not need refinement
                // else{
                //     // Generate the particle inside the node
                    
                //     // Exception for the nodes in baseGrid (Still put in the idle list)
                //     if (IDflag.at(_ID)) continue;

                //     // *** Here is the problem (Problem occured here)
                //     if (baseGrid.nodeMap.count(_Node->headNodeID) == 0) continue; //  WHY THIS IS NEEDED ???

                //     // Create aliasing to HEAD node
                //     Node *&_headNode = baseGrid.nodeMap.at(_Node->headNodeID);

                //     // Create a particle neighbor candidate ID list
                //     std::vector<Node*> nodeNghIDList;
                //     std::vector<int> parNghIDList;
                        
                //     // Find all LEAF neighbor at the HEAD node
                //     baseGrid.findNghAll(nodeNghIDList, _headNode);
                //     for (auto &_nghNode : nodeNghIDList){
                //         // // Push all particle into the neighbor particle ID list
                //         // for (auto &_nghParID : _nghNode->parList){
                //         //     parNghIDList.push_back(_nghParID);
                //         // }
                //         if (parNodeMap.count(_nghNode->nodeID) != 0){     // Make sure the lagrangian particle is inside the node
                //             for (auto &_nghParID : parNodeMap.at(_nghNode->nodeID)){
                //                 parNghIDList.push_back(_nghParID);
                //             }
                //         }
                //     }

                //     // Generate the particle
                //     double _headSize = _headNode->length / this->parNum;
                //     this->generateParticle(_Node, _headSize);

                //     // Find the neighbor
                //     if (nghConv == 0)
                //     this->findParticleNeighbor(_Node->parList, parNghIDList, currPar);
                // }
            } /* LOOP [3] -> Check at each neighbor */

            // Release the queue container at the current level
            for (int &_ID : IDqueue[currLvl]){
                IDflag.erase(_ID);          // Release the ID flag
            }
            IDqueue[currLvl].clear();       // Release the ID list

        } /* LOOP [2] -> Check at each level */

        // // DEBUG LINE
        // std::string name = "AdaptationNLD" + std::to_string(iter++);
        // baseGrid.saveLeafGrid(tempGrid, name);
        // save_particle(*this->newPar, name);
        // // iter++;
    
    } /* LOOP [1] -> Loop through all queue */

    // _time = clock() - _time;
    // printf("<+> Number of refined node              : %8d\n", refinedNode);
    // printf("<-> Node refinement computational time [%f s]\n\n", (double)_time/CLOCKS_PER_SEC);
    
    // Display calculation time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<+> Adapt 3: NLD Check     : %f s\n", _time);

    // // [DEBUG] Check the temporary grid
    // baseGrid.saveLeafGrid(tempGrid,"AdaptationFinNLD");
    // baseGrid.saveGrid(tempGrid,"AdaptationFinNLD");

    // // ========================================
    // // ------------ THIRD DEBUGGER ------------
    // // ========================================
    // std::cout << "\n[LOG 3] PRINT NUMBER OF NODE!\n";
    // std::cout << FONT_GRAY;
    // std::cout << " >>> Refinement  : " << this->refineList.size() << "\n";
    // std::cout << " >>> IDLE        : " << this->idleList.size() << "\n";
    // std::cout << " >>> Compression : " << this->compressList.size() << "\n\n";
    // std::cout << FONT_RESET;
        

    // PROCEDURE 4!
    // ************
    // Compress all nodes in compression list: each node is called a HEAD node:
    /* Do the compression
        > Check sibling inside the map
           -> If not : Put all sibling to idle map
           -> If yes proceed
        > Then check all sibling have to compression
           -> If not : Put all sibling to idle map
           -> If yes proceed
        > Check the NLD (Neighbor level different) criteria
           -> Find all in level current node parent neighbor
           -> Check whether the leaf of each neighbor is inside idle, refine or compress list
               -> if inside IDLE do basic NLD
               -> if inside refine check the leaf at refine
               -> if inside compression (continue to next ngh List)
    */
    
    // Compression node ID leveled queue container
    std::vector<std::vector<int>> toCompressIDList(this->maxLevel+1);
    // Assign the node ID into the queue container
    for (auto &[_ID,_flag] : this->compressList){
        Node *&_node = baseGrid.nodeMap.at(_ID);
        toCompressIDList[_node->level].push_back(_ID);
    }

    // START PROCEDURE 4:
    int compressNode = 0;        // Compression counter for console log
    
    // MESSAGE_LOG << "Check the compression criteria\n";
    iter = 0;
	
    // Time counter
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif
    
    // At this stage, all compression node in temporary GRID is not existed yet
    // *[1] Perform node compression
    // [LOOP 1] -> Evaluate the node from the finest level
    for (int currLev = this->maxLevel; currLev > 0; currLev--){
        // Loop is not including the rood node (ROOT node cannot be compressed)
        // std::cout << FONT_GREEN << " >>> Print the current level : " << currLev << "\n" << FONT_RESET;

        // [LOOP 2] -> Iterate through all node in the compression ID in the list
        for (auto &_ID : toCompressIDList[currLev]){
            // At initial toCompressIDList contain all node to be compressed
            // But further, the element may include
            //    [1] parent of compressed node and
            //    [2] non compressed node
            // A flow control is needed to separate the node (to be executed)

            // Aliasing the current node
            Node *_node = baseGrid.nodeMap.at(_ID);     // Based on the based Grid
            
            // Check whether the current node need to be compressed
            if (!(_node->needCompression)) continue;    // See the bottom algorithm to know why <!>

            // * Compression checking procedure, the compression will continue if
            // [1] All siblings need to be compressed (inside the compression List)
            // [2] The parent of compressed node meet the NLD criteria

            
            // **[CHECK 1] -> Check sibling
            // Find the sibling node list
            std::vector<Node*> siblingList;
            baseGrid.findSibling(siblingList, _node);

            // Check whether all sibling are to be compressed
            bool proceed = true;
            for (auto _sibNode : siblingList){
                // Check whether the sibling is need to be compressed
                if (!(_sibNode->needCompression)){
                    // The current sibling is not need to compress
                    proceed = false;
                    break;
                }
            }
            
            // Not all sibling need to compress
            if (!proceed){
                // Shut down all compression flag in this sibling
                for (auto _sibNode : siblingList){
                    // Shut down all the compression flag
                    _sibNode->needCompression = false;

                    // *For Particle generation
                    // Do modification only for the sibling inside the compression list and a LEAF NODE
                    if (this->compressList.count(_sibNode->nodeID) != 0){
                        // Move the node from compression list into the idle list
                        if (_sibNode->isLeaf){      // [!] I'm not sure why need to check for leaf <?>
                            // Close the compression flag at this sibling node
                            this->compressList.erase(_sibNode->nodeID);
                            // Put into the idle list
                            this->idleList.insert({_sibNode->nodeID,true});
                        }
                    }
                }
                continue;   // Proceed to the next iteration
            }


            // **[CHECK 2] -> Check the NLD criteria on the parent
            // At this stage, all sibling are known inside the compression list
            
            // Initialize the parent node parameter
            int _parNodeID = baseGrid.findParent(_node);        // Parent ID
            Node *_parNode = baseGrid.nodeMap.at(_parNodeID);   // Parent node
            int _parLevel = _parNode->level;                    // Parent level

            // Find the neighbor node list of the parent node
            std::vector<Node*> allLeafNeighborNode;
            _parNode->isLeaf = true;
            baseGrid.findNghAll(allLeafNeighborNode, _parNode);
            _parNode->isLeaf = false;
            
            // MISSING POINT: The neighbor node maybe located on the refinement list !

            // <!> A custom addition code line to assign the active sign
            // [!] Addition of the Shut down all the compression flag
            for (auto _sibNode : siblingList){
                if (_sibNode->isActive){
                    _parNode->isActive = true;
                    break;
                }
            }

            // Initialize a compression execution flag
            bool needCompress = true;

            // [LOOP for CHECK 2] -> Iterate through all LEAF neighbor node of parent node
            for (auto &_nghNode : allLeafNeighborNode){
                // This loop only need to check the level difference between parent and all neighbors
                // are no more than MAX_LEVEL_DIFF

                // No need to check the current node
                if ((_nghNode->nodeID) == _parNodeID) continue;

                // A potential of the neighbor ID
                // * inside COMPRESSION -> Check if already compressed
                // * inside REFINEMENT  -> Check the maximum Level
                // * inside IDLE        -> Direct evaluation

                // NOTE: I supposed the following condition should using an if, else if, else structure
                // **The neighbor lies on the compression list
                if (this->compressList.count(_nghNode->nodeID) != 0){
                    // Aliasing to the neighbor node
                    const Node *_tempNode;

                    // Determine on which grid the node is located
                    if(this->tempGrid.nodeMap.count(_nghNode->nodeID) != 0){
                        // If the node is existed in the temporary node, It must be the previously compressed node
                        _tempNode = this->tempGrid.nodeMap.at(_nghNode->nodeID);
                        
                        // A recursive to the parent node if not a LEAF node
                        do {
                            // Find the parent ID
                            int __parID = baseGrid.findParent(_tempNode);
                            _tempNode = this->tempGrid.nodeMap.at(__parID);  // Change the temporary node to its parent

                            // Skip the loop if fullfil the NDL criteria
                            if ((_tempNode->level) <= (_parLevel + this->NghlevelDiff)) break;
                        } while (!(_tempNode->isLeaf));

                    }
                    else{
                        // This node have not undergo any compression yet
                        _tempNode = baseGrid.nodeMap.at(_nghNode->nodeID);
                    }

                    // Direct check to a LEAF node, only check if the compress flag still on
                    if ((_tempNode->level) > (_parLevel + this->NghlevelDiff)){
                        // This node doesn't fulfill the NLD criteria
                        needCompress = false;
                        break;
                    }

                    if(needCompress) continue;
                    else break;
                }

                // **The neighbor lies on the refine list
                if (this->refineList.count(_nghNode->nodeID) != 0){
                    // Proceed to child and child until the end
                    std::vector<int> _evalIDList;   // A queue list of all child
                    std::vector<int> _tempChdList;  // Temporary child ID container

                    // Push the current neighbor as the first evaluated node
                    _evalIDList.push_back(_nghNode->nodeID);
                    
                    // Iterate all node to be evaluated
                    for (size_t i = 0; i < _evalIDList.size(); i++){
                        // Aliasing the node
                        Node *_evalNode = this->tempGrid.nodeMap.at(_evalIDList[i]);
                        int _evalLvl = _evalNode->level;
                        
                        // Check the NLD criteria
                        if (!_evalNode->isLeaf){
                            // Check the NLD criteria of the child
                            if (_evalLvl + 1 > _parLevel + this->NghlevelDiff){
                                // This node doesn't fulfill the NLD criteria
                                needCompress = false;
                                break;
                            }
                            // Still fullfil the NLD criteria -> navigate to child
                            _tempChdList = baseGrid.findChild(_evalNode);
                            for (auto &_chdID : _tempChdList){
                                _evalIDList.push_back(_chdID);
                            }
                        }
                    }
                    
                    if (needCompress) continue;
                    else break;
                }
                
                // **The neighbor lies on the idle list
                if ((_nghNode->level) > (_parLevel + this->NghlevelDiff)){
                    // This node doesn't fulfill the NLD criteria
                    needCompress = false;
                    break;
                }
            }

            // *Execute compression
            if (needCompress){
                // * MODIFICATION OF ALL SIBLING NODE
                for (size_t i = 0; i < siblingList.size(); i++){
                    // Create duplicate node of sibling to the temporary Grid (Only for the non existing node)
                    if (this->tempGrid.nodeMap.count(siblingList[i]->nodeID) == 0){
                        // Duplicate the sibling node into the temporary grid
                        Node *_dupNode = new Node(siblingList[i]);
                        _dupNode->isLeaf = false;

                        // Push the duplicate node into the temporary grid
                        this->tempGrid.nodeMap.insert({_dupNode->nodeID, _dupNode});
                    }else{
                        this->tempGrid.nodeMap.at(siblingList[i]->nodeID)->isLeaf = false;
                    }
                    
                    // Turn off the compression flag for all sibling list
                    siblingList[i]->needCompression = false;

                    // Turn of the compression flag at the compression container of the sibling node
                    this->compressList.at(siblingList[i]->nodeID) = false;  // <?> NEEDED FOR PARTICLE GENERATION <?>
                }
                
                // Set the target resolution level of the parent node
                for (auto _sibNode : siblingList) _parNode->tarResLv = 
                std::max<int>(_parNode->tarResLv, _sibNode->tarResLv);
                
                // Set up the parent node compression flag (both node in temporary and base grid)
                _parNode->needCompression = _parNode->tarResLv < _parNode->level ? true : false;
                
                // * GENERATE THE PARENT NODE (new node after compression)
                // Create duplicate node of parent node
                Node *_dupParNode = new Node(_parNode);
                _dupParNode->isLeaf = true;
                compressNode++;      // Node is compressed add into account

                // * PUSH NODE INTO CONTAINER
                this->tempGrid.nodeMap.insert({_parNodeID, _dupParNode});   // The grid Node
                this->compressList.insert({_parNodeID, true});              // The compression flag list <?> NEEDED FOR PARTICLE GENERATION <?>
                toCompressIDList[_parLevel].push_back(_parNodeID);          // The queue list
            }
            
            // *Not execute compression
            else{
                // std::cout << FONT_DARK_BLUE << "[DEBUG] SOMETHING WEIRD HAPPENED ON THE NODE " <<_ID << FONT_RESET << "\n";
                // Shut down all compression flag on all this sibling
                for (auto _sibNode : siblingList){
                    // Shut down all the compression flag
                    _sibNode->needCompression = false;

                    // *For particle generation
                    // Check if the sibling node is a leaf
                    if (_sibNode->isLeaf){
                        // Close the compression flag of sibling
                        this->compressList.erase(_sibNode->nodeID);
                        // Put the sibling into the idle list
                        this->idleList.insert({_sibNode->nodeID, true});
                    }
                }
            }
        }

        // // DEBUG LINE saving the data
        // std::string name = "COMPRESSION" + std::to_string(iter++);
        // baseGrid.saveLeafGrid(tempGrid, name);
    }

    // Display calculation time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<+> Adapt 4: Compression 1 : %f s\n", _time);

    // Compress the current NODE
    /* Last Note:
        > Create the node [Only the node], in the temporary grid list (Put the headnodeID as the first child ID)
        > Put node into the compressed container
        > HOW to hanlDe multicompression (?) : The old ID is multiple ?? 
            -> Create a node division funciton like
            -> This function put all generated particle inside all leaf of current cmopression node.
            -> List all leaf node in this compressed Node -> To check the neighbor of all particle inside each leaf
            -> NOTE THAT all leaf is inside the base grid
        -> 
    */

    // MESSAGE_LOG << "Generate particle inside compressed node\n";
    // Time counter
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif
	
    // // *[2] Generate particel at each compressed node
    // for (auto &[_currID, _currFlag] : this->compressList){
    //     // Only compress a true compressed List
    //     if (_currFlag){    
    //         // Create an alias to this node
    //         Node *_currNode = this->tempGrid.nodeMap.at(_currID);

    //         // Generate the particle on the current node
    //         double _headSize = _currNode->length / this->parNum;
    //         this->generateParticle(_currNode, _headSize);

    //         // The node ID container of each child node from the current node
    //         std::vector<int> _leafNodeIDList;               // Leaf node list that curr on this node
    //         std::vector<int> _queueNodeIDList = {_currID};  // A queue container for particle division
            
    //         // *Promptly divide (inherit) the particle into the very LEAF node
    //         for (size_t i = 0; i < _queueNodeIDList.size(); i++){
    //             // Internal variable
    //             int _evalID = _queueNodeIDList[i];
                
    //             // Check LEAF node of the current node (check toward the base grid)
    //             if (baseGrid.nodeMap.at(_evalID)->isLeaf){
    //                 _leafNodeIDList.push_back(_evalID);
    //                 continue;
    //             }
                
    //             // Find the children of the current node
    //             Node *_evalNode = this->tempGrid.nodeMap.at(_evalID);
    //             std::vector<int> _currChdIDList = baseGrid.findChild(_evalNode);

    //             // Put all children into the queue container
    //             for (int &_chdID : _currChdIDList){
    //                 _queueNodeIDList.push_back(_chdID);
    //             }

    //             // Divide the particle into new child node
    //             this->tempGrid.divideParticle(_evalNode, _newPar, _currChdIDList, 0);
    //         }
            

    //         // *Evaluate the neighbor of all LEAF node
    //         for (int &_evalID : _leafNodeIDList){
    //             // Create an alias
    //             Node *_evalNode = this->tempGrid.nodeMap.at(_evalID);
    //             double _tarSize = _evalNode->length / this->parNum;
                
    //             // Update the particle size
    //             for (int &_parID : _evalNode->parList) _newPar.s[_parID] = _tarSize;

    //             // Find all neighbor id based on the current before refinement ID
    //             std::vector<Node*> nodeNghIDList;
    //             std::vector<int> parNghIDList;  // List of all neighbor particle candidate
    //             baseGrid.findNghAll(nodeNghIDList, baseGrid.nodeMap.at(_evalID)); // Find all neighbor at the old node location
    //             for (auto &_nghNode : nodeNghIDList){
    //                 // // Push all particle inside the current node into the ngh par ID list (Inlcude the old node too)
    //                 // for (auto &_nghParID : _nghNode->parList){
    //                 //     parNghIDList.push_back(_nghParID);
    //                 // }
    //                 if (parNodeMap.count(_nghNode->nodeID) != 0){     // Make sure the lagrangian particle is inside the node
    //                     for (auto &_nghParID : parNodeMap.at(_nghNode->nodeID)){
    //                         parNghIDList.push_back(_nghParID);
    //                     }
    //                 }
    //             }
                
    //             // Evaluate the neighbor
    //             if (nghConv == 0)
    //             this->findParticleNeighbor(_evalNode->parList, parNghIDList, currPar);
    //         }
    //     }
    // }

    // Display calculation time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<+> Adapt 5: Compression 2 : %f s\n", _time);

    // // DEBUG LINE
    // baseGrid.saveLeafGrid(tempGrid, "AdaptationCompressed");
    // save_particle(*this->newPar, "AdaptationCompressed");

    // _time = clock() - _time;
    // printf("<+> Number of compressed node           : %8d\n", compressNode);
    // printf("<-> Node compress computational time   [%f s]\n\n", (double)_time/CLOCKS_PER_SEC);

    // // ========================================
    // // ------------ FOURTH DEBUGGER -----------
    // // ========================================
    // std::cout << "\n[LOG 4] PRINT NUMBER OF NODE!\n";
    // std::cout << FONT_GRAY;
    // std::cout << " >>> Refinement  : " << this->refineList.size() << "\n";
    // std::cout << " >>> IDLE        : " << this->idleList.size() << "\n";
    // std::cout << " >>> Compression : " << this->compressList.size() << "\n\n";
    // std::cout << FONT_RESET;


    // No need adaptation if there is no refinement and compression
    if (this->refineList.empty() && this->compressList.empty()){
        std::cout << FONT_MAROON << "[FINAL LOG] "
                  << FONT_RESET  << "Be carefull, adaptation is cancelled!!!\n";
        // Prevent memory leak
        delete this->newPar;
        return;
    }

    
    // PROCEDURE 5!
    // ************
    // Update the particle and the grid node data
    // [1] Insert the particle data from idle List
    // [2] Combine the node data from temporary list

    // MESSAGE_LOG << "[DATA 1] Merge the particle data!!\n";
    // Time counter
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif

    // // *[1] Take all particle inside idle list into the new partilce
    // for (auto &[_ID,_flag] : this->idleList){
    //     Node *&_node = baseGrid.nodeMap.at(_ID);
    //     // Get the list of particle inside the node
    //     std::vector<int> tempParList = _node->parList;

    //     // Also update the particle data of the current node
    //     _node->parList.clear();
        
    //     // Transfer all particle data from the particle list inside the node
    //     for (int &_particleID : tempParList){
    //         _newPar.x.push_back(tempPar.x[_particleID]);
    //         _newPar.y.push_back(tempPar.y[_particleID]);
    //         if (DIM > 2) _newPar.z.push_back(tempPar.z[_particleID]);
    //         _newPar.s.push_back(tempPar.s[_particleID]);
    //         _newPar.isActive.push_back(tempPar.isActive[_particleID]);
    //         _newPar.neighbor.push_back(std::vector<int>());
    //         // if (nghConv == 0){
    //         //     _newPar.neighbor.push_back(tempPar.neighbor[_particleID]);
    //         //     if (!Pars::flag_ngh_include_self){
    //         //         // Put the current ID into the neighbor list
    //         //         _newPar.neighbor[_newPar.num].push_back(_particleID);
    //         //     }
    //         // }
    //         _newPar.level.push_back(tempPar.level[_particleID]);
    //         _newPar.nodeID.push_back(tempPar.nodeID[_particleID]);
            
    //         // Update the particle list in the node
    //         _node->parList.push_back(_newPar.num);

    //         _newPar.num++;
    //     }

    //     // printf("[1] Done particle data assignment on node %d\n", _ID);

    //     // Find the current node neighbor node consisted of leaf node only
    //     // Create a particle neighbor candidate ID list
    //     // std::vector<Node*> tempNodeNghIDList;
    //     std::vector<Node*> nodeNghIDList;
    //     std::vector<int> parNghIDList;
    //     // std::unordered_map<int, bool> parentNodeCnt;    // Count list of the parent node into list

    //     // Find all LEAF neighbor at the HEAD node
    //     baseGrid.findNghAll(nodeNghIDList, _node);

    //     // printf("[2] Done neighbor calculation on node %d\n", _ID);
        
    //     // Check whether the node is included in the temporary grid or not
    //     // for (auto &_nghNode : tempNodeNghIDList){
    //     //     if (this->tempGrid.nodeMap.count(_nghNode->nodeID) != 0){
    //     //         // Put an alias to the current node
    //     //         Node *&_node = this->tempGrid.nodeMap.at(_nghNode->nodeID);

    //     //         // Here supposed to find the parent or child node
    //     //         if (this->refineList.count(_nghNode->nodeID) != 0){
    //     //             // A refinary list
    //     //             //  -> Find child until meet the leaf node
    //     //             std::vector<int> chdIDList;
    //     //             chdIDList = baseGrid.findChild(_node);
                    
    //     //             // A recursive container
    //     //             std::vector<Node*> nodeQueue;
    //     //             for (auto& chdID : chdIDList) nodeQueue.push_back(this->tempGrid.nodeMap.at(chdID));

    //     //             // Do recursively until find the leaf
    //     //             for (auto &_evalNode : nodeQueue){
    //     //                 // Evaluate the current node, whether a leaf or not
    //     //                 // Check whether the parent is a leaf or not
    //     //                 if (_evalNode->isLeaf){
    //     //                     nodeNghIDList.push_back(_evalNode);
    //     //                 }else{
    //     //                     // Evaluate the child again
    //     //                     chdIDList = baseGrid.findChild(_evalNode);
    //     //                     for (auto& chdID : chdIDList) nodeQueue.push_back(this->tempGrid.nodeMap.at(chdID));
    //     //                 }
    //     //             }
    //     //         }else{
    //     //             // A compression list
    //     //             //  -> Find parent until meet the leaf
    //     //             std::vector<Node*> nodeQueue = {_node};
                    
    //     //             // Do recursively until find the leaf
    //     //             for (auto &_evalNode : nodeQueue){
    //     //                 // Evaluate parent node ID
    //     //                 int nghParID = baseGrid.findParent(_node);

    //     //                 // Check whether the parent already evaluated or not
    //     //                 if (parentNodeCnt.count(nghParID) == 0){   // Still not existed in the list
    //     //                     // Update the container
    //     //                     parentNodeCnt.insert({nghParID, true});

    //     //                     // Alias to the parent node
    //     //                     Node*& parNode = this->tempGrid.nodeMap.at(nghParID);

    //     //                     // Check whether the parent is a leaf or not
    //     //                     if (parNode->isLeaf){
    //     //                         nodeNghIDList.push_back(parNode);
    //     //                     }else{
    //     //                         nodeQueue.push_back(parNode);
    //     //                     }
    //     //                 }
    //     //             }
    //     //         }
    //     //     }else{
    //     //         // Located in idle list
    //     //         nodeNghIDList.push_back(_nghNode);
    //     //     }
    //     // }

    //     // Push all particle into the neighbor particle ID list
    //     for (auto &_nghNode : nodeNghIDList){
    //         // // Push all particle into the neighbor particle ID list
    //         // for (auto &_nghParticleID : _nghNode->parList){
    //         //     parNghIDList.push_back(_nghParticleID);
    //         // }
    //         if (parNodeMap.count(_nghNode->nodeID) != 0){     // Make sure the lagrangian particle is inside the node
    //             for (auto &_nghParID : parNodeMap.at(_nghNode->nodeID)){
    //                 parNghIDList.push_back(_nghParID);
    //             }
    //         }
    //     }

    //     // printf("[3] Done evaluation on the node neighbor node %d\n", _ID);

    //     // Evaluate the neighbor
    //     if (nghConv == 0)
    //     this->findParticleNeighbor(_node->parList, parNghIDList, currPar);

    //     // printf("[FINAL] Done neighbor calculation on node %d\n", _ID);
    // }

    // printf("<+> DONE PARTICLE GENERATION\n");
    
    // DEBUG LINE DEBUG LINE DEBUG LINE DEBUG LINE
    // ===========================================
    // // Saving Address of the node
    // std::ofstream writter;
    // writter.open("address_data_tmp.csv", std::ofstream::out);
    // writter << "ID,Add\n";
    // for (auto &[key,val] : this->tempGrid.nodeMap){
    //     writter << key 
    //             << "," << val
    //             << "\n";
    // }
    // writter.close();

    // writter.open("address_data_base.csv", std::ofstream::out);
    // writter << "ID,Add\n";
    // for (auto &[key,val] : baseGrid.nodeMap){
    //     writter << key 
    //             << "," << val
    //             << "\n";
    // }
    // writter.close();

    // MESSAGE_LOG << "[DATA 2] Update the grid node data!!\n";
    // *[2] Update the node list of base grid from temporary grid
    for (auto &[_ID,_flag] : this->tempGrid.nodeMap){
        // Check if the current Node is existed in the base grid map
        if (baseGrid.nodeMap.count(_ID) != 0){
            // Existed at the base grid node we need to delete the node on baseGrid

            // Check whether the node is in the refine or in the compression head
            if (this->refineList.count(_ID) != 0 || this->compressList.at(_ID)){
                delete baseGrid.nodeMap.at(_ID);    // Free the memory
                baseGrid.nodeMap.erase(_ID);        // Delete from the container
                // Delete the current node, replace with the one in the temporary grid
            }else{
                // Delete the node in both Grid Node
                delete baseGrid.nodeMap.at(_ID);        // Free the memory
                baseGrid.nodeMap.erase(_ID);    // Delete from the container
                delete this->tempGrid.nodeMap.at(_ID);  // Free the memory
                continue;
            }
        }
        
        // Push the current node into the map
        Node *&_node = this->tempGrid.nodeMap.at(_ID);
        baseGrid.nodeMap.insert({_ID, _node});
    }

    // // DEBUG LINE
    // baseGrid.saveGrid(tempGrid, "AdaptationDone");
    // save_particle(*this->newPar, "AdaptationDone");

    // Remove all ID from tempGrid
    // for (auto &[key,val] : this->tempGrid.nodeMap) delete val;
    this->tempGrid.nodeMap.clear();     // Here we share the node to the base grid Node

    // Display calculation time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    
    // // *END OF PROCEDURE 5*
    // baseGrid.saveLeafGrid(baseGrid, "AdaptationFinal");
    // baseGrid.saveGrid(baseGrid, "AdaptationFinal");
    
    // // DEBUG: Show the number of particle neighbor
    // for (int i = 0; i < this->newPar->num; i++){
    //     std::cout << "Particle " << i << ": neighbor :" << this->newPar->neighbor[i].size() << "\n";
    //     // for (int j = 0; j < this->newPar->neighbor[i].size(); j++){
    //     //     std::cout << 
    //     // }
    //     // if (i == 1000) break;
    // }

    // // DEBUG: Save the data of neighbor
    // // Save the data neighbor data of particle 
    // std::ofstream writeData;   
    // // std::vector<int> dataID = {338032, 342838, 338038, 309182, 300361};
    // // std::vector<int> dataID = {342805, 345224, 346830, 322846, 312378};
    // std::vector<int> dataID = {0, 10, 300, 245, 125};
    // for (auto ID : dataID){
    //     writeData.open("output/CollectNgh_"+std::to_string(ID)+".csv");
    //     writeData << "x,y,z,s\n";
    //     // Collect the grid Data first
    //     writeData << this->newPar->x[ID]
    //               << "," << this->newPar->y[ID]
    //               << "," << this->newPar->z[ID]
    //               << "," << this->newPar->s[ID]
    //               << "\n";
    //     for (int nghID : this->newPar->neighbor[ID]){
    //     writeData << currPar.x[nghID]
    //               << "," << currPar.y[nghID]
    //               << "," << currPar.z[nghID]
    //               << "," << currPar.s[nghID]
    //               << "\n";
    //     }
    //     writeData.close();
    // }


    // Calculate the total leaf node
    int leafCnt = 0;
    for (auto& [_ID, _node] : baseGrid.nodeMap){
        if (_node->isLeaf) leafCnt++;
    }
    MESSAGE_LOG << "Final total nodes      : " << baseGrid.nodeMap.size() << "\n";
    MESSAGE_LOG << "Final total leaf nodes : " << leafCnt << "\n";
    MESSAGE_LOG << "Final total particle   : " << this->newPar->num << "\n";
    std::cout << FONT_GREEN << "[FINAL LOG] Adaptation executed successfully!!!" << FONT_RESET << "\n" ;
    
    // The adaptation is performed
    return;
}
