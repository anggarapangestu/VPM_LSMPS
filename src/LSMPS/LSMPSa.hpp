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

class LSMPSa
{
private:
    // Internal variables
    const double R_fac = 3.1;      // effective radius ratio
    static const int MAT_SIZE = 5; // size of matrix

    std::vector<double> _ddx;    // d{}/dx
    std::vector<double> _ddy;    // d{}/dy
    std::vector<double> _d2d2x;  // d^2{}/d^2x
    std::vector<double> _d2dxdy; // d^2{}/dxdy
    std::vector<double> _d2d2y;  // d^2{}/d^2y

    // LSMPS A Utilities
    std::vector<double> get_p(const double &x, const double &y, const double &s);
    double weight_function(const double &rij, const double &Rij);

    // Direct LSMPS A calculation
    void calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                         const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist);
    // LSMPS A with collocation calculation
    void calculate_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &sSource, const std::vector<double> &fSource,
                         const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                         std::vector<std::vector<int>> &neighborlist);
public:
    // Direct LSMPS A
    void set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                   const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist);
    // LSMPS A with collocation (Diffusion)
    void set_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &sSource, const std::vector<double> &fSource,
                   const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                   std::vector<std::vector<int>> &neighborlist);

    // Spatial differential getter
    std::vector<double> get_ddx();
    std::vector<double> get_ddy();
    std::vector<double> get_d2d2x();
    std::vector<double> get_d2dxdy();
    std::vector<double> get_d2d2y();
};
#endif