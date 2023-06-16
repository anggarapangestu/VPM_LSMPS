#ifndef INCLUDED_PRESSURE_POISSON
#define INCLUDED_PRESSURE_POISSON

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

#ifndef INCLUDED_LSMPSa
#include "../LSMPS/LSMPSa.hpp"
#endif

#ifndef INCLUDED_LSMPSb
#include "../LSMPS/LSMPSb.hpp"
#endif

#ifndef FMM_2D_PACKAGE
#include "../FMM2D/fmm2D.hpp"
#endif

class pressure_poisson{
private:
    // -- creating instance
	treeCell treeData;  // the basis of tree data of poisson solver
    int counter = 0;
    int counter1 = 0;
    void save_pressure_source(const std::vector<double>&x, const std::vector<double>&b, const std::vector<double>&src, const std::vector<bool>&mark);
public:
    void get_pressure(Particle& par);
};

#endif
