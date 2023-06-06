#ifndef INCLUDED_LSMPSa
#define INCLUDED_LSMPSa
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <omp.h>          // pragma omp ...
#include "../Eigen/Dense" // linear algebra library
#include "../Eigen/Sparse"
using namespace Eigen;    // for linear algebra operation

//// template <typename T>
class LSMPSa
{
private:
//// #define MAT_SIZE 5 // size of matrix
    const double R_fac = 3.1;      // effective radius ratio
    static const int MAT_SIZE = 5; // size of matrix

    std::vector<double> _ddx;    // d{}/dx
    std::vector<double> _ddy;    // d{}/dy
    std::vector<double> _d2d2x;  // d^2{}/d^2x
    std::vector<double> _d2dxdy; // d^2{}/dxdy
    std::vector<double> _d2d2y;  // d^2{}/d^2y

    std::vector<double> get_p(const double &x, const double &y, const double &s);
    double weight_function(const double &rij, const double &Rij);

    void calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                         const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist);
    void calculate_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &sSource, const std::vector<double> &fSource,
                         const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                         std::vector<std::vector<int>> &neighborlist);
    //void calculate_LSMPS_Laplace(std::vector<double> &u, std::vector<double> &v, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
    //                     std::vector<double> &f, std::vector<std::vector<int>> &neighborlist, std::vector<int> isboundary);

    /*void calculate_ETA_LSMPS(std::vector<std::vector<double>> &LSMPS_EtaDx2, std::vector<std::vector<double>> &LSMPS_EtaDy2, std::vector<double> x, std::vector<double> y, 
                        std::vector<std::vector<int>> neighbour, std::vector<double> e);*/
public:
    // LSMPS(/* args */);
    // ~LSMPS();

    void set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                   const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist);
    void set_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &sSource, const std::vector<double> &fSource,
                   const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                   std::vector<std::vector<int>> &neighborlist);
    //void set_LSMPS_Laplace(std::vector<double> &u,  std::vector<double> &v, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
    //                    std::vector<double> &f, std::vector<std::vector<int>> &neighborlist, std::vector<int> boundary);

    std::vector<double> get_ddx();
    std::vector<double> get_ddy();
    std::vector<double> get_d2d2x();
    std::vector<double> get_d2dxdy();
    std::vector<double> get_d2d2y();
};
#endif