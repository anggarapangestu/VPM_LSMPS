#ifndef PARTICLE_GRID_NODE_ADAPTATION
#define PARTICLE_GRID_NODE_ADAPTATION

#include "../../Utils.hpp"
#include "../grid_block/gridNode.hpp"

/**
 *  @brief  Adaptation class consisted of method to evaluate and  
 *  perform the particle distribution adaptation inside the
 *  simulation domain.
 * 
 *  @headerfile remeshing.hpp
*/
class GridNodeAdapt
{
private:
    // Class member variable
    std::vector<double> maxValue;        // The maximum value to evaluate target resolution level (TRL) of a node
    const double winTol = std::log(Pars::adapt_tol);
    double maxLevel;        // The maximum level at this "GridNode" evaluation
    int NghlevelDiff;       // NLD criteria of maximum neighbor difference level
    std::unordered_map<int,bool> refineList;    // All node ID list to be refine
    std::unordered_map<int,bool> compressList;  // All node ID list to be compression
    std::unordered_map<int,bool> idleList;      // All node ID list that not doing anything

    GridNode tempGrid;      // Temporary GridNode to hold the node that undergo refinement or compression
    Particle* newPar;       // New particle distribution after adaptation

    // Private member for particle generation
    int parNum;         // The number of particle in a node in one dimension
    int totParNum;      // The total number of particle in a node
    int div[DIM];       // Iterator divisor for particle generation

    // INITIALIZATION METHOD
    void targetLevel(GridNode &grid, const Particle &par, const std::vector<std::vector<double>*>& errPredProp);
    void targetLevelv2(GridNode &grid, const Particle &par, const std::vector<std::vector<double>*>& errPredProp,
                        std::unordered_map<int,std::vector<int>> &GridCont);
    void assignPar2Grid(const GridNode &grid, const Particle &par, std::unordered_map<int,std::vector<int>> &GridCont);
    int targetLevelCalc(Node *_node, int parID, const std::vector<std::vector<double>*>& errPredProp);
    void set_adaptation_flag(GridNode &grid);

    // PROCEDURE 1&2 METHOD {Actually nothing happened here}
    int assignParticle(Node *currNode, const Particle &templatePar, const std::vector<int> &chdIDlist, int type);

    void generateParticle(Node *currNode, double parSize);
    // void generateParticleActive(Node *currNode, double parSize, bool activeSign);
    void findParticleNeighbor(const std::vector<int> &parIDList, const std::vector<int> &parNghIDList, const Particle &nghSrcPar);

public:
    // Adaptation process

    bool get_adaptation(GridNode &baseGrid, const Particle &tempPar, 
                        const Particle &currPar, int nghConv,
                        std::vector<std::vector<double>*>& errPredProp,
                        std::unordered_map<int, std::vector<int>> &parNodeMap);

    bool get_adaptation(GridNode &baseGrid, const Particle &tempPar, 
                        const Particle &currPar,
                        std::vector<std::vector<double>*>& errPredProp);
    bool get_adaptation(GridNode &baseGrid, const Particle &tempPar, 
                        const Particle &currPar, int nghConv);
    void take_particle_pointer(Particle *&_tarPar);

    void get_adaptation_LSMPS(GridNode &baseGrid, const Particle &currPar, 
                              std::vector<std::vector<double>*>& errPredProp);
    
    // Default Constructor
    GridNodeAdapt(): maxValue(0.0){
        // Initialize the value of predefined member
        this->NghlevelDiff = Pars::ngh_diff_level;
        this->maxLevel = Pars::max_level;
        
        // Reserve the dynamic container
        refineList.clear();
        compressList.clear();
        idleList.clear();
    }
    
    // Default Deconstructor
    ~GridNodeAdapt(){};
};


#endif