#ifndef INCLUDED_NEIGHBOR_GRID_NODE
#define INCLUDED_NEIGHBOR_GRID_NODE

#include "gridNode.hpp"

class GridNodeNgh
{
public:
    void find_neighbor(std::vector<std::vector<int>> &_nghIDList,
                       const GridNode &_baseGrid,
                       const Particle &_evalPar);
    void find_inter_neighbor(std::vector<std::vector<int>> &_nghIDList,
                             const Particle &_evalPar,
                             const std::unordered_map<int, std::vector<int>> &_evalParNodeMap,
                             const GridNode &_baseGrid,
                             const Particle &_sourcePar);

    void eval_inter_ngh_gridNode(std::vector<std::vector<int>> &_nghIDList,
                                 const Particle &_trgPar,
                                 const Particle &_srcPar,
                                 const GridNode &_trgGridNode,
                                 const std::unordered_map<int, std::vector<int>> &_srcParNodeMap);

    void assign_par2node(const GridNode &baseGrid, 
                         std::unordered_map<int, std::vector<int>> &parNodeMap, 
                         Particle &evalPar);
};

#endif