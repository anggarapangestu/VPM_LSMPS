#ifndef INCLUDED_LSMPSb
#define INCLUDED_LSMPSb
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <omp.h>          // pragma omp ...
#include "../Eigen/Dense" // linear algebra library

using namespace Eigen;    // for linear algebra operation

//// template <typename T>
class LSMPSb
{
private:
//// #define MAT_SIZE 6 // size of matrix
    const double R_fac = 2.1;       // effective radius ratio awalnya 3.1
    static const int MAT_SIZE = 6;  // size of matrix

    std::vector<double> _d00;    // d^
    std::vector<double> _ddx;    // d{}/dx
    std::vector<double> _ddy;    // d{}/dy
    std::vector<double> _d2d2x;  // d^2{}/d^2x
    std::vector<double> _d2dxdy; // d^2{}/dxdy
    std::vector<double> _d2d2y;  // d^2{}/d^2y
    std::vector<double> _d3d3x;  // d^3{}/d^3x

    std::vector<double> get_p(const double &x, const double &y, const double &s);
    double weight_function(const double &rij, const double &Rij);

    void calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s, const std::vector<double> &f,
                         const std::vector<double> &xi, const std::vector<double> &yi, const std::vector<double> &si, const std::vector<double> &fi,
                         std::vector<std::vector<int>> &neighborlistInter, std::vector<std::vector<int>> &neighborlistIntra);

public:
    void set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s, const std::vector<double> &f,
                   const std::vector<double> &xi, const std::vector<double> &yi, const std::vector<double> &si, const std::vector<double> &fi,
                   std::vector<std::vector<int>> &neighborlistInter, std::vector<std::vector<int>> &neighborlistIntra);

    std::vector<double> get_d00();
    std::vector<double> get_ddx();
    std::vector<double> get_ddy();
    std::vector<double> get_d2d2x();
    std::vector<double> get_d2dxdy();
    std::vector<double> get_d2d2y();
};
#endif