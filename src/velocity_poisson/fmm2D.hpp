#ifndef FMM_2D_PACKAGE
#define FMM_2D_PACKAGE

#include "../../global.hpp"

#ifndef INCLUDED_TREE_CELL
#include "treeCell.hpp"
#endif

class fmm2D{
private:
    // Internal variable
    int parNum;         // Number of particle
    int expOrd;         // Maximum expansion order
    int cellNum;        // Total number of cell 
    int maxLevel;       // Highest tree level

    // The FMM result data of each particle
    std::vector<double> parPotential;   // The result data potential
    std::vector<double> parField_x;     // The result data field in x direction
    std::vector<double> parField_y;     // The result data field in y direction
    
    // FMM Source Expansion Value (The particle source sum; indexing -> [cellID, expansion order])
    std::vector<std::vector<std::complex<double>>> ak;  // Local particle source sum, origin at current cell center
    std::vector<std::vector<std::complex<double>>> bk;  // Total particle source sum, origin at current cell center

    // Internal method
    int binComb(int n, int k);  // Function to calculate the binomial combinatoric
    void potDirSum(int _cellID, treeCell& cellData, std::vector<int>& _nghCell, 
                   std::vector<std::vector<double>>& pos, std::vector<double>& src);    // Function to calculate the direct sum potential
    void fieldDirSum(int _cellID, treeCell& cellData, std::vector<int>& _nghCell, 
                     std::vector<std::vector<double>>& pos, std::vector<double>& src);  // Function to calculate the direct sum field

public:
    // The FMM calculation for potential
    void setupFMM(treeCell& cellData, std::vector<std::vector<double>>& pos, 
    std::vector<bool>& activeMark, std::vector<double>& src);

    // The FMM calculation for potential
    void calcPotential(treeCell& cellData, std::vector<std::vector<double>>& pos, 
    std::vector<bool>& activeMark, std::vector<double>& src);
    
    // The FMM calculation for field
    void calcField(treeCell& cellData, std::vector<std::vector<double>>& pos, 
    std::vector<bool>& activeMark, std::vector<double>& src);

    // The FMM get function from the result of calculation
    void get_Potential(std::vector<double>& phi);
    void get_Field(std::vector<double>& Ex, std::vector<double>& Ey);
};

#endif