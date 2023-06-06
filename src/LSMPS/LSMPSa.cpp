#include "LSMPSa.hpp"

#pragma region public methods
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void LSMPSa::set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                       const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist)
{
    // clear local variables
    this->_ddx.clear();
    this->_ddy.clear();
    this->_d2d2x.clear();
    this->_d2dxdy.clear();
    this->_d2d2y.clear();

    // resize local variabels
    int nparticle = x.size();
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2d2y.resize(nparticle);

    // perform calculation
    this->calculate_LSMPS(x, y, s, f, neighborlist);
    
    {   // TODO: perform testing
        // double a = -9.97681;
        // double b = -5.28855;
        // double c = 2.96304;
        // double d = -8.51253;
        // double e = -4.59517;
        // ////double f = -2.79946;
        // std::ofstream ofs10, ofs01, ofs11, ofs20, ofs02;
        // ofs10.open("output/lsmps_10.csv");
        // ofs01.open("output/lsmps_01.csv");
        // ofs11.open("output/lsmps_11.csv");
        // ofs20.open("output/lsmps_20.csv");
        // ofs02.open("output/lsmps_02.csv");
        // for (size_t i = 0; i < nparticle; i++)
        // {
        //     double _ddxAnalytic = 2 * a * x[i] + c * y[i] + d;
        //     double _ddyAnalytic = 2 * b * y[i] + c * x[i] + e;
        //     double _d2d2xAnalytic = 2 * a;
        //     double _d2dxdyAnalytic = c;
        //     double _d2d2yAnalytic = 2 * b;

        //     double _ratio10 = std::abs(_ddx[i]) / std::abs(_ddxAnalytic);
        //     double _ratio01 = std::abs(_ddy[i]) / std::abs(_ddyAnalytic);
        //     double _ratio11 = std::abs(_d2dxdy[i]) / std::abs(_d2dxdyAnalytic);
        //     double _ratio20 = std::abs(_d2d2x[i]) / std::abs(_d2d2xAnalytic);
        //     double _ratio02 = std::abs(_d2d2y[i]) / std::abs(_d2d2yAnalytic);

        //     ofs10 << x[i] << "," << y[i] << "," << _ddxAnalytic << "," << _ddx[i] << "," << _ratio10 << "\n";
        //     ofs01 << x[i] << "," << y[i] << "," << _ddyAnalytic << "," << _ddy[i] << "," << _ratio01 << "\n";
        //     ofs11 << x[i] << "," << y[i] << "," << _d2dxdyAnalytic << "," << _d2dxdy[i] << "," << _ratio11 << "\n";
        //     ofs20 << x[i] << "," << y[i] << "," << _d2d2xAnalytic << "," << _d2d2x[i] << "," << _ratio20 << "\n";
        //     ofs02 << x[i] << "," << y[i] << "," << _d2d2yAnalytic << "," << _d2d2y[i] << "," << _ratio02 << "\n";
        // }
        // ofs10.close();
        // ofs01.close();
        // ofs11.close();
        // ofs20.close();
        // ofs02.close();
    }
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void LSMPSa::set_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &sSource, const std::vector<double> &fSource,
                       const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                       std::vector<std::vector<int>> &neighborlist)
{
    // clear local variables
    this->_ddx.clear();
    this->_ddy.clear();
    this->_d2d2x.clear();
    this->_d2dxdy.clear();
    this->_d2d2y.clear();

    // resize local variabels
    int nparticle = xSource.size();
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2d2y.resize(nparticle);

    // perform calculation
    this->calculate_LSMPS(xSource, ySource, sSource, fSource, xCollocation, yCollocation, sCollocation, fCollocation, neighborlist);
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_ddx()
{
    return this->_ddx;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_ddy()
{
    return this->_ddy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2d2x()
{
    return this->_d2d2x;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2dxdy()
{
    return this->_d2dxdy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2d2y()
{
    return this->_d2d2y;
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
#pragma endregion

#pragma region private methods
// ===========================================================================
// ===========================================================================
void LSMPSa::calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                             const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist)
{
    int nparticle = x.size();
    
    // Evaluate the LSMPS A of all particle
    for (size_t i = 0; i < nparticle; i++)
    {
        // Initialize the LSMPS matrices
        Eigen::MatrixXd Hi = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::MatrixXd Mi = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::VectorXd bi = Eigen::VectorXd::Zero(MAT_SIZE);
        
        // Initialize H_rs Matrix
        Hi(0, 0) = std::pow(s[i], -1);
        Hi(1, 1) = std::pow(s[i], -1);
        Hi(2, 2) = std::pow(s[i], -2) * 2;
        Hi(3, 3) = std::pow(s[i], -2);
        Hi(4, 4) = std::pow(s[i], -2) * 2;

        // Evaluate at each neighbor particle
        std::vector<int> _neighborIndex = neighborlist[i];
        for (size_t j = 0; j < _neighborIndex.size(); j++)
        {
            // Note that neighbor particle and evaluated particle inside the same variable
            int idxi = i;                   // Particle ID
            int idxj = _neighborIndex[j];   // Neighbor Particle ID
            
            // Particle position
            double _xi = x[idxi];
            double _yi = y[idxi];
            
            // Neighbor particle position
            double _xj = x[idxj];
            double _yj = y[idxj];

            // Coordinate distance between Particle and Neighbor Particle
            double _xij = _xj - _xi;
            double _yij = _yj - _yi;
            
            // LSMPS A interpolated variable value different
            double _fij = f[idxj] - f[idxi];    // Neighbor particle - evaluated particle

            // Resultant distance(_rij) and effective radius(_Rij)
            double _rij = std::sqrt(std::pow(_xij, 2) + std::pow(_yij, 2)); //distance between particles. 
            double _Ri = this->R_fac * s[idxi];   // Effective radius of Current particle
            // double _Rj = this->R_fac * s[idxj];   // Effective radius of Neighbor particle
            // double _Rij = (_Ri + _Rj) * 0.5;      // Effective radius of Average particle 

            std::vector<double> _p1 = this->get_p(_xij, _yij, s[idxi]);
            std::vector<double> _p2 = this->get_p(_xij, _yij, s[idxi]);
            double _wij = this->weight_function(_rij, _Ri) * std::pow(s[idxj]/s[idxi],2);

            // Calculatin moment matrix M and b
            for (size_t k1 = 0; k1 < MAT_SIZE; k1++)
            {
                for (size_t k2 = 0; k2 < MAT_SIZE; k2++)
                {
                    // generate tensor product between p
                    Mi(k1, k2) = Mi(k1, k2) + (_wij * _p1[k1] * _p2[k2]);
                }
                // generate moment matrix
                bi(k1) = bi(k1) + (_wij * _p1[k1] * _fij);
            }
        }

        // Solve Least Square
        Eigen::VectorXd MiInv_Bi = Mi.bdcSvd(ComputeThinU | ComputeThinV).solve(bi); // (MAT_SIZE x 1)
        Eigen::VectorXd Dx = Hi * MiInv_Bi;                                          // (MAT_SIZE x 1)

        // Assign to private variables
        this->_ddx[i] = Dx[0];
        this->_ddy[i] = Dx[1];
        this->_d2d2x[i] = Dx[2];
        this->_d2dxdy[i] = Dx[3];
        this->_d2d2y[i] = Dx[4];
    }
}
// ===========================================================================
// ===========================================================================
void LSMPSa::calculate_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &sSource, const std::vector<double> &fSource,
                             const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                             std::vector<std::vector<int>> &neighborlist)
{
    /*
    LSMPS A collocation Note :
    > There are two relative position in this operation
      * The particle target interpolation (refer as SOURCE)
      * The particle interpolation data (refer as COLLOCATION)
    > Please note that SOURCE and COLLOCATION just a reference name
      * Both SOURCE and COLLOCATION is particle variable
      * Both SOURCE and COLLOCATION come from the same particle variable
    
    SOURCE VARIABLE - interpolation target
    vec<dbl> xSource, vec<dbl> ySource  : Particle position list
    vec<dbl> sSource                    : Particle size list
    vec<dbl> fSource                    : Interpolated variable data

    COLLOCATION VARIABLE - interpolation data source
    vec<dbl> xCollocation, vec<dbl> yCollocation   : Particle position list
    vec<dbl> sCollocation                          : Particle size list
    vec<dbl> fCollocation                          : Interpolated variable data
    vec<vec<int>> neighborlist                     : Grid neighbor ID list consisted of particle ID
    */

    int nparticle = xSource.size();

    // Evaluate the LSMPS A of all particle
    for (size_t i = 0; i < nparticle; i++)
    {
        // Initialize the LSMPS matrices
        Eigen::MatrixXd Hi = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::MatrixXd Mi = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::VectorXd bi = Eigen::VectorXd::Zero(MAT_SIZE);
        
        // Initialize H_rs Matrix
        Hi(0, 0) = std::pow(sSource[i], -1);
        Hi(1, 1) = std::pow(sSource[i], -1);
        Hi(2, 2) = std::pow(sSource[i], -2) * 2;
        Hi(3, 3) = std::pow(sSource[i], -2);
        Hi(4, 4) = std::pow(sSource[i], -2) * 2;

        // Evaluate at each neighbor particle
        std::vector<int> _neighborIndex = neighborlist[i];
        for (size_t j = 0; j < _neighborIndex.size(); j++)
        {
            // Note that neighbor particle and evaluated particle inside the same variable
            int idxi = i;                      // Particle ID
            int idxj = _neighborIndex[j];      // Neighbor Particle ID

            // SOURCE particle position (Evaluated particle)
            double _xi = xSource[idxi];
            double _yi = ySource[idxi];

            // COLLOCATION particle position (Neighbor particle)
            double _xj = xCollocation[idxj];
            double _yj = yCollocation[idxj];

            // Coordinate distance between SOURCE and COLLOCATION Particle
            double _xij = _xj - _xi;
            double _yij = _yj - _yi;
            
            // LSMPS A interpolated variable value different
            double _fij = fCollocation[idxj] - fSource[idxi];
            
            // Resultant distance(_rij) and effective radius(_Rij)
            double _rij = std::sqrt(std::pow(_xij, 2) + std::pow(_yij, 2));
            double _Ri = this->R_fac * sSource[idxi];         // Effective radius of Current particle
            // double _Rj = this->R_fac * sCollocation[idxj];    // Effective radius of Neighbor particle
            // double _Rij = (_Ri + _Rj) * 0.5;                  // Effective radius of Average particle

            std::vector<double> _p1 = this->get_p(_xij, _yij, sSource[idxi]);
            std::vector<double> _p2 = this->get_p(_xij, _yij, sSource[idxi]);
            double _wij = this->weight_function(_rij, _Ri) * std::pow(sCollocation[idxj]/sSource[idxi],2);

            // Calculation of moment matrix M and b
            for (size_t k1 = 0; k1 < MAT_SIZE; k1++)
            {
                for (size_t k2 = 0; k2 < MAT_SIZE; k2++)
                {
                    // Generate tensor product between p
                    Mi(k1, k2) = Mi(k1, k2) + (_wij * _p1[k1] * _p2[k2]);
                }
                // Generate moment matrix
                bi(k1) = bi(k1) + (_wij * _p1[k1] * _fij);
            }
        }

        // Solve Least Square
        Eigen::VectorXd MiInv_Bi = Mi.bdcSvd(ComputeThinU | ComputeThinV).solve(bi); // (MAT_SIZE x 1)
        Eigen::VectorXd Dx = Hi * MiInv_Bi;                                          // (MAT_SIZE x 1)

        // Assign to private variables
        this->_ddx[i] = Dx[0];
        this->_ddy[i] = Dx[1];
        this->_d2d2x[i] = Dx[2];
        this->_d2dxdy[i] = Dx[3];
        this->_d2d2y[i] = Dx[4];
    }
}
// ===========================================================================
// ===========================================================================
std::vector<double> LSMPSa::get_p(const double &xij, const double &yij, const double &si)
{
    std::vector<double> _p(MAT_SIZE);

    double _xij = xij / si;
    double _yij = yij / si;

    _p[0] = _xij;
    _p[1] = _yij;
    _p[2] = _xij * _xij;
    _p[3] = _xij * _yij;
    _p[4] = _yij * _yij;

    return _p;
}
// ===========================================================================
// ===========================================================================
double LSMPSa::weight_function(const double &rij, const double &Rij)
{
    double _wij;
    if (rij <= Rij)
    {
        _wij = std::pow(1 - (rij / Rij), 2);
    }
    else
    {
        _wij = 0.0e0;
    }

    return _wij;
}
// ===========================================================================
// ===========================================================================
/*void LSMPSa::set_LSMPS_Laplace(std::vector<double> &u, std::vector<double> &v, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                        std::vector<double> &f, std::vector<std::vector<int>> &neighborlist, std::vector<int> boundary)
{
    // clear local variab_ddx.clear();
    _ddy.clear();
    _d2d2x.clear();
    _d2dxdy.clear();
    _d2d2y.clear();

    // resize local variabels
    int nparticle = x.size();
    //_ddx.resize(nparticle);
    //_ddy.resize(nparticle);
    //_d2d2x.resize(nparticle);
    //_d2dxdy.resize(nparticle);
    //_d2d2y.resize(nparticle);

    // perform calculation
    calculate_LSMPS_Laplace(u, v, x, y, s, f, neighborlist, boundary);
}
*/
#pragma endregion
/*
void LSMPSa::calculate_LSMPS_Laplace(std::vector<double> &u, std::vector<double> &v, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                             std::vector<double> &f, std::vector<std::vector<int>> &neighborlist, std::vector<int> isboundary)
{
    int nparticle = x.size();
    std::vector<std::vector<double>> EtaDxu;
    std::vector<std::vector<double>> EtaDyu;
    std::vector<std::vector<double>> EtaDxv;
    std::vector<std::vector<double>> EtaDyv;
    std::vector<Eigen::MatrixXd> HrsMinvData(nparticle);
    std::vector<std::vector<Eigen::MatrixXd>> bdata(nparticle);


    for (size_t i = 0; i < nparticle; i++)
    {
        std::vector<int> _neighborIndex = neighborlist[i];

        Eigen::MatrixXd Hi = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::MatrixXd M = Eigen::MatrixXd::Zero(MAT_SIZE,MAT_SIZE);
        std::vector<Eigen::MatrixXd> btemp(_neighborIndex.size());
        std::vector<double> weight(_neighborIndex.size());
        Eigen::MatrixXd P(MAT_SIZE,1);

        Hi(0, 0) = std::pow(s[i], -1);
        Hi(1, 1) = std::pow(s[i], -1);
        Hi(2, 2) = std::pow(s[i], -2) * 2;
        Hi(3, 3) = std::pow(s[i], -2);
        Hi(4, 4) = std::pow(s[i], -2) * 2;

        
        for (size_t j = 0; j < _neighborIndex.size(); j++)
        {
            int idxi = i;
            int idxj = _neighborIndex[j];
            double _xi = x[idxi];
            double _yi = y[idxi];
            double _xj = x[idxj];
            double _yj = y[idxj];
            double _xij = _xj - _xi;
            double _yij = _yj - _yi;

            double _rij = std::sqrt(std::pow(_xij, 2) + std::pow(_yij, 2)); 
            double _Ri = R_e * s[idxi];
            double _Rj = R_e * s[idxi];
            double _Rij = (_Ri + _Rj) * 0.5;
            weight[j] = weight_function(_rij, _Rij);
        }

        for (size_t j = 0; j < _neighborIndex.size(); j++)
        {
            int idxi = i;
            int idxj = _neighborIndex[j];
            double _xi = x[idxi];
            double _yi = y[idxi];
            double _xj = x[idxj];
            double _yj = y[idxj];
            double _xij = _xj - _xi;
            double _yij = _yj - _yi;
            
            _xij = _xij / s[idxi];
            _yij = _yij / s[idxj];

            P(0,0) = _xij;
            P(1,0) = _yij;
            P(2,0) = _xij*_xij;
            P(3,0) = _xij*_yij;
            P(4,0) = _yij*_yij;

            M = M + weight[j]*(P*P.transpose());
            btemp[j] = weight[j]*P;
        }

        Eigen::MatrixXd Minv = M.inverse();
        Eigen::MatrixXd HrsMinv = Hi*Minv;

        HrsMinvData[i] = HrsMinv;
        bdata[i] = btemp;

        for (size_t j = 0; j < _neighborIndex.size(); j++)
        {
            int jdx = _neighborIndex[j];

            MatrixXd Eta_temp = (HrsMinvData[i]*bdata[i][j]);

            EtaDxu[i].push_back(Eta_temp(1,0)); 
            EtaDyu[i].push_back(Eta_temp(2,0));
        }
    }

    //clear vector
    HrsMinvData.resize(0);
    bdata.resize(0);
    
    //Openmp Threads
    setNbThreads(4);

    // Source Term
    Eigen::VectorXd b = Eigen::VectorXd::Zero(2*nparticle);
    for (int i = 0; i < nparticle; i++)
    {
        b(i) = 0;
        b((2*i)+1) = f[i];
    }

    // Matrix Term


    //for nabla cross u = w
    typedef Eigen::Triplet<double> T2;
    std::vector<T2> tripletlist2;

    //for nabla dot u = 0
    typedef Eigen::Triplet<double> T1;
    std::vector<T1> tripletlist1;

    for(size_t i = 0; i < nparticle; i++)
    {
        // Nabla.U
        int row;
        int col;
        double data1, data2; //data1 for ux, data2 for uy
        double temp1, temp2; //temp1 for ux, temp2 for uy
        
        row = i;
        temp1 = 0.0;
        temp2 = 0.0;

        for (int j = 0; j < neighborlist[i].size(); j++)
        {
            col = neighborlist[i][j];
            data1 = EtaDxu[i][j];
            data2 = EtaDyu[i][j];
            tripletlist1.push_back(T1(row,2*col,data1));
            tripletlist1.push_back(T1(row,2*col+1,data2));
            temp1 = temp1 - data1;
            temp2 = temp2 - data2;
        }
        col = i;
        data1 = temp1;
        data2 = temp2;
        tripletlist1.push_back(T1(row,2*col,data1));
        tripletlist1.push_back(T1(row,2*col+1,data2));

        //nabla cross u 
        if (isboundary[i] == 1)
        {
            row = i;
            col = i;
            data = 1.0;
            tripletlist.push_back(T(row,col,data));
        }
        else
        {
            row = i;
            temp = 0.0;
            for (int j = 0; j < neighborlist[i].size(); j++)
            {
                col = neighborlist[i][j];
                data = (EtaDxu[i][j] + EtaDyu[i][j]);
                tripletlist.push_back(T(row,col,data));
                temp = temp - data;
            }
            col = i;
            data = temp;
            tripletlist.push_back(T(row,col,data));
        }
    }

    Eigen::SparseMatrix<double> A (2*nparticle, 2*nparticle);
    A.setFromTriplets(tripletlist.begin(),tripletlist.end());

    A.makeCompressed();
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    Eigen::VectorXd u_eig = solver.solve(b);

}*/

