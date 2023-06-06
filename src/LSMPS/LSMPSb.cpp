#include "LSMPSb.hpp"

#pragma region public methods
// ===========================================================================
// ***************************************************************************
void LSMPSb::set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s, const std::vector<double> &f,
                       const std::vector<double> &xi, const std::vector<double> &yi, const std::vector<double> &si, const std::vector<double> &fi,
                       std::vector<std::vector<int>> &neighborlistInter, std::vector<std::vector<int>> &neighborlistIntra)
{
    // clear local variables
    this->_d00.clear();
    this->_ddx.clear();
    this->_ddy.clear();
    this->_d2d2x.clear();
    this->_d2dxdy.clear();
    this->_d2d2y.clear();
    this->_d3d3x.clear();

    // resize local variabels
    int nparticle = x.size();
    this->_d00.resize(nparticle);
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2d2y.resize(nparticle);
    this->_d3d3x.resize(nparticle);

    // perform calculation
    this->calculate_LSMPS(x, y, s, f, xi, yi, si, fi, neighborlistInter, neighborlistIntra);

// TODO: perform testing
#pragma region testing
    // double _a = -9.97681;
    // double _b = -5.28855;
    // double _c = 2.96304;
    // double _d = -8.51253;
    // double _e = -4.59517;
    // double _f = -2.79946;
    // std::ofstream ofs00, ofs10, ofs01, ofs11, ofs20, ofs02;
    // ofs00.open("output/lsmpsb_00.dat");
    // ofs10.open("output/lsmpsb_10.dat");
    // ofs01.open("output/lsmpsb_01.dat");
    // ofs11.open("output/lsmpsb_11.dat");
    // ofs20.open("output/lsmpsb_20.dat");
    // ofs02.open("output/lsmpsb_02.dat");
    // for (size_t i = 0; i < nparticle; i++)
    // {
    //     double _analytic = _a * std::pow(x[i], 2) + _b * std::pow(y[i], 2) + _c * x[i] * y[i] + _d * x[i] + _e * y[i] + _f;
    //     double _ddxAnalytic = 2 * _a * x[i] + _c * y[i] + _d;
    //     double _ddyAnalytic = 2 * _b * y[i] + _c * x[i] + _e;
    //     double _d2d2xAnalytic = 2 * _a;
    //     double _d2dxdyAnalytic = _c;
    //     double _d2d2yAnalytic = 2 * _b;

    //     double _ratio00 = std::abs(_d00[i]) / std::abs(_analytic);
    //     double _ratio10 = std::abs(_ddx[i]) / std::abs(_ddxAnalytic);
    //     double _ratio01 = std::abs(_ddy[i]) / std::abs(_ddyAnalytic);
    //     double _ratio11 = std::abs(_d2dxdy[i]) / std::abs(_d2dxdyAnalytic);
    //     double _ratio20 = std::abs(_d2d2x[i]) / std::abs(_d2d2xAnalytic);
    //     double _ratio02 = std::abs(_d2d2y[i]) / std::abs(_d2d2yAnalytic);

    //     ofs00 << x[i] << " " << y[i] << " " << _analytic << " " << _d00[i] << " " << _ratio00 << "\n";
    //     ofs10 << x[i] << " " << y[i] << " " << _ddxAnalytic << " " << _ddx[i] << " " << _ratio10 << "\n";
    //     ofs01 << x[i] << " " << y[i] << " " << _ddyAnalytic << " " << _ddy[i] << " " << _ratio01 << "\n";
    //     ofs11 << x[i] << " " << y[i] << " " << _d2dxdyAnalytic << " " << _d2dxdy[i] << " " << _ratio11 << "\n";
    //     ofs20 << x[i] << " " << y[i] << " " << _d2d2xAnalytic << " " << _d2d2x[i] << " " << _ratio20 << "\n";
    //     ofs02 << x[i] << " " << y[i] << " " << _d2d2yAnalytic << " " << _d2d2y[i] << " " << _ratio02 << "\n";
    // }
    // ofs00.close();
    // ofs10.close();
    // ofs01.close();
    // ofs11.close();
    // ofs20.close();
    // ofs02.close();
#pragma endregion
}
// ***************************************************************************
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d00()
{
    return this->_d00;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_ddx()
{
    return this->_ddx;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_ddy()
{
    return this->_ddy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d2d2x()
{
    return this->_d2d2x;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d2dxdy()
{
    return this->_d2dxdy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d2d2y()
{
    return this->_d2d2y;
}
// ---------------------------------------------------------------------------
// ***************************************************************************
#pragma endregion

#pragma region private methods
// ===========================================================================
// ***************************************************************************
void LSMPSb::calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s, const std::vector<double> &f,
                             const std::vector<double> &xi, const std::vector<double> &yi, const std::vector<double> &si, const std::vector<double> &fi,
                             std::vector<std::vector<int>> &neighborlistInter, std::vector<std::vector<int>> &neighborlistIntra)
{
    /*
    LSMPS B Note :
    > There are two relative position in this operation
      * The target interpolated position (refer as GRID)
      * The position of interpolation data (refer as PARTICLE)
    > Please note that GRID and PARTICLE just a reference name
      * Both GRID and PARTICLE is particle variable
    
    GRID VARIABLE - interpolation target
    vec<dbl> x, vec<dbl> y  : Particle position list
    vec<dbl> s              : Particle size list
    vec<dbl> f              : Interpolated variable data
    
    PARTICLE VARIABLE - interpolation data source
    vec<dbl> xi, vec<dbl> yi    : Particle position list
    vec<dbl> si                 : Particle size list
    vec<dbl> fi                 : Interpolated variable data
    vec<vec<int>> neighborlistInter : Grid neighbor ID list consisted of particle ID
    vec<vec<int>> neighborlistIntra : Grid neighbor ID list consisted of grid ID
    */
    
    int nparticle = x.size();
    int nClear = 0;
    int nLow = 0;
    std::vector<double> averageDiameter(nparticle, 0.0e0);

    // Evaluate the LSMPS B of all GRID particle (NOT the PARTICLE variable !)
    for (size_t i = 0; i < nparticle; i++)
    {
        // Initialize the LSMPS matrices
        Eigen::MatrixXd Hbar = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::MatrixXd Mbar = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::VectorXd bbar = Eigen::VectorXd::Zero(MAT_SIZE);

        // Evaluate the neighbor of the evaluated GRID
        std::vector<int> _neighborInterIndex = neighborlistInter[i];
        std::vector<int> _neighborIntraIndex = neighborlistIntra[i];

        if (_neighborInterIndex.empty())
        {
            // If there is no nearby particle from the grid, the grid differential value set to zero
            this->_d00[i] = 0;
            this->_ddx[i] = 0;
            this->_ddy[i] = 0;
            this->_d2d2x[i] = 0;
            this->_d2dxdy[i] = 0;
            this->_d2d2y[i] = 0;
            this->_d3d3x[i] = 0;
        }
        else
        {
            // // [OLD CODE]
            // // TODO: calculate average diameter
            // for (size_t j = 0; j < _neighborInterIndex.size(); j++)
            // {
            //     int idx = _neighborInterIndex[j];
            //     averageDiameter[i] += si[idx];
            // }
            
            // // Why using this as r_s <?>
            // averageDiameter[i] = averageDiameter[i] / (double)_neighborInterIndex.size();
            
            // Use the GRID size
            averageDiameter[i] = s[i];
            
            // Initialize H_rs Matrix
            Hbar(0, 0) = 1;
            Hbar(1, 1) = std::pow(averageDiameter[i], -1);
            Hbar(2, 2) = std::pow(averageDiameter[i], -1);
            Hbar(3, 3) = std::pow(averageDiameter[i], -2) * 2;
            Hbar(4, 4) = std::pow(averageDiameter[i], -2);
            Hbar(5, 5) = std::pow(averageDiameter[i], -2) * 2;
            //Hbar(6, 6) = std::pow(averageDiameter[i], -3) * 6;

            // TODO: perform Least Square MPS: Evaluate each neighbor particle
            for (size_t j = 0; j < _neighborInterIndex.size(); j++)
            {
                int idxi = i;                       // GRID ID
                int idxj = _neighborInterIndex[j];  // PARTICLE ID as neighbor
                
                // GRID position
                double _xi = x[idxi]; 
                double _yi = y[idxi];
                
                // Neigbor: PARTICLE position
                double _xj = xi[idxj];
                double _yj = yi[idxj];
                
                // Coordinate distance between GRID and PARTICLE
                double _xij = _xj - _xi;       // x distance between GRID[i] and neighbor: PARTICLE[j]
                double _yij = _yj - _yi;       // y distance between GRID[i] and neighbor: PARTICLE[j]
                
                // Interpolated variable of neighbor: PARTICLE
                double _fij = fi[idxj] - 0.0e0; 
                
                // Resultant distance(_rij) and effective radius(_Rij)
                double _rij = std::sqrt(std::pow(_xij, 2) + std::pow(_yij, 2)); 
                double _Ri = this->R_fac * averageDiameter[idxi];   // Effective radius of GRID particle
                // double _Rj = this->R_fac * si[idxj];                // Effective radius of PARTICLE particle
                // double _Rij = (_Ri + _Rj) * 0.5;                    // Effective radius of Average particle
                
                // Calculate p and weight
                std::vector<double> _p1 = this->get_p(_xij, _yij, averageDiameter[i]);  
                std::vector<double> _p2 = this->get_p(_xij, _yij, averageDiameter[i]);  
                double _wij = this->weight_function(_rij, _Ri) * std::pow(si[idxj]/averageDiameter[idxi],2);

                // Calculate the moment matrix M and b
                for (size_t k1 = 0; k1 < MAT_SIZE; k1++)
                {
                    for (size_t k2 = 0; k2 < MAT_SIZE; k2++)
                    {
                        // generate tensor product between p
                        Mbar(k1, k2) = Mbar(k1, k2) + (_wij * _p1[k1] * _p2[k2]);
                    }
                    // generate moment matrix
                    bbar(k1) = bbar(k1) + (_wij * _p1[k1] * _fij);
                }
            }
            
            // TODO: to increase robustness <- considering intra-particle neighbor
            if (std::abs(Mbar.determinant()) < 1.0e-6) // ! if tend to be NON-invertible matrix
            {
                // FOR DEBUGGING
                nClear --;
                nLow ++;
                // END FOR DEBUGGING
                for (size_t j = 0; j < _neighborIntraIndex.size(); j++)
                {
                    int idxi = i;                        // GRID ID
                    int idxj = _neighborIntraIndex[j];   // GRID ID as neighbor
                    
                    // GRID position
                    double _xi = x[idxi];
                    double _yi = y[idxi];
                    
                    // Neigbor: GRID position
                    double _xj = x[idxj];
                    double _yj = y[idxj];

                    // Coordinate distance between GRID and GRID:Neighbor
                    double _xij = _xj - _xi;       // x distance between GRID[i] and neighbor: GRID[j]
                    double _yij = _yj - _yi;       // y distance between GRID[i] and neighbor: GRID[j]
                    
                    // Interpolated variable of neighbor: GRID
                    double _fij = 0.0e0;           // ! no interaction of the function, value set to 0.0, then no need to do LSMPS<?>

                    // Resultant distance(_rij) and effective radius(_Rij)
                    double _rij = std::sqrt(std::pow(_xij, 2) + std::pow(_yij, 2));
                    double _Ri = this->R_fac * averageDiameter[idxi];   // Effective radius of PARTICLE particle
                    // double _Rj = this->R_fac * s[idxj];                 // Effective radius of PARTICLE particle
                    // double _Rij = (_Ri + _Rj) * 0.5;                    // Effective radius of PARTICLE particle

                    // Calculate p and weight
                    std::vector<double> _p1 = this->get_p(_xij, _yij, averageDiameter[i]);
                    std::vector<double> _p2 = this->get_p(_xij, _yij, averageDiameter[i]);
                    double _wij = this->weight_function(_rij, _Ri) * std::pow(s[idxj]/s[idxi],2);

                    // Calculate the moment matrix M and b
                    for (size_t k1 = 0; k1 < MAT_SIZE; k1++)
                    {
                        for (size_t k2 = 0; k2 < MAT_SIZE; k2++)
                        {
                            // generate tensor product between p
                            Mbar(k1, k2) = Mbar(k1, k2) + (_wij * _p1[k1] * _p2[k2]);
                        }
                        // generate moment matrix
                        bbar(k1) = bbar(k1) + (_wij * _p1[k1] * _fij);
                    }
                }
            }

            // Solve Least Square Matrix Operation
            // --- Method 1 ---
            // LLT<MatrixXd> llt;
            // llt.compute(Mbar);
            // Eigen::VectorXd MbarInv_Bbar = llt.solve(bbar);
            // --- Method 2 ---
            Eigen::VectorXd MbarInv_Bbar = Mbar.fullPivHouseholderQr().solve(bbar);
            // std::cout << "[DEBUG] Matrix Calculated\n";
            // --- Method 3 ---
            // Eigen::VectorXd MbarInv_Bbar = Mbar.bdcSvd(ComputeThinU | ComputeThinV).solve(bbar); // (MAT_SIZE x 1)
            
            // Final result
            Eigen::VectorXd Dx = Hbar * MbarInv_Bbar; // (MAT_SIZE x 1)

            // Assign the LSMPS B result
            this->_d00[i] = Dx[0];
            this->_ddx[i] = Dx[1];
            this->_ddy[i] = Dx[2];
            this->_d2d2x[i] = Dx[3];
            this->_d2dxdy[i] = Dx[4];
            this->_d2d2y[i] = Dx[5];
            //this->_d3d3x[i] = Dx[6];
            nClear++;
        }
    }
    printf("<+> Number of clear interpolation       : %8d\n", nClear);
    printf("<+> Number of low determinant           : %8d\n", nLow);
}
// ===========================================================================
// ***************************************************************************
std::vector<double> LSMPSb::get_p(const double &xij, const double &yij, const double &si)
{
    std::vector<double> _p(MAT_SIZE);

    double _xij = xij / si;
    double _yij = yij / si;

    _p[0] = 1;
    _p[1] = _xij;
    _p[2] = _yij;
    _p[3] = _xij * _xij;
    _p[4] = _xij * _yij;
    _p[5] = _yij * _yij;
    //_p[6] = _xij * _xij * _xij;

    return _p;
}
// ===========================================================================
// ***************************************************************************
double LSMPSb::weight_function(const double &rij, const double &Rij)
{
    double _wij;
    if (rij <= Rij)
    {
        // Quadratic Weight Function
        _wij = std::pow(1 - (rij / Rij), 2);
        
        // // Singular Weight Function
        // _wij = std::pow((rij / Rij) - 1, 2);
    }
    else
    {
        _wij = 0.0e0;
    }

    return _wij;
}
// ===========================================================================
// ***************************************************************************
#pragma endregion