#ifndef INCLUDED_TREE_CELL
#define INCLUDED_TREE_CELL

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

#include <iostream>
#include <sstream>
#include <fstream>

class treeCell{
private:
    // TREE BASIC PARAMETERs
    double baseLen;     // The length of the root cell
    double expansion;   // The tree basis domain expansion factor
    int max_level;      // The maximum level of the tree
    int min_level;      // The minimum level of leaf cell
    int max_par;        // The maximum number of source particle in cell
    int basis;          // The number of cell division
    int cellNum;        // The total number of the cell in the tree
    
    // Transformation map : in determining the neighboring cell
    std::vector<std::vector<int>> cellIdx;              // The corresponding transformation matrix index
    std::vector<std::vector<std::vector<int>>> T_map2D; // The transformation matrix cell IDs for 2D
    std::vector<std::vector<std::vector<std::vector<int>>>> T_map3D; // The transformation matrix cell IDs for 3D
    
    void createTmap2D();  // Creating the Tree map for 2D
    void createTmap3D();  // Creating the Tree map for 3D

public:
    // CELL DATA VARIABLE : with ID in sequence of the parent - child
    std::vector<std::vector<double>> cellPos;   // The center position of the cell
    std::vector<bool> leafCellMark;             // The marking for leaf cell or not
    std::vector<bool> outsideCell;              // The marking whether the cell is existed but outside the domain
    std::vector<int> level;                     // The cell level hierarchy of the cell
    std::vector<int> parNum;                    // The number of all particel inside the cell
    std::vector<int> srcNum;                    // The number of source particel inside the cell
    std::vector<std::vector<int>> parIDList;    // List of all particle ID in the current cell

    // Additional for FMM calculation
    std::vector<int> leafCellList;              // List of all leaf cell

    // External method
    int findParentID(int currID);
    std::vector<int> findChildID(int currID);

    // Basic External method
    void initializeTree(std::vector<std::vector<double>>& parPos);
    void createTree(std::vector<std::vector<double>>& parPos, std::vector<bool> & activeMark);

    void updateCell(std::vector<std::vector<double>>& parPos, std::vector<bool> & activeMark);

    /* Find the neighbor list
       4 type neighbor cell
       * Type 1 : Touched neighbor, leaf cell
       * Type 2 : Well separated neighbor (current parent touched neighbor parent cell), at the same level to current cell
       * Type 3 : Well separated neighbor (current cell touched parent), leaf cell
       * Type 4 : Well separated neighbor (current parent touched neighbor parent cell), at low level to current cell
    */

    void intList(int currID, std::vector<int>& ID_list_1, std::vector<int>& ID_list_3);
    void extList(int currID, std::vector<int>& ID_list_2, std::vector<int>& ID_list_4);
    void findLevel(const std::vector<int>& ID_list, std::vector<int>& level);
    
    // Save data
    void saveTreeCell();
    void saveTreeCell(std::string name);
};

#endif