#ifndef INCLUDED_UTILS		// Already include global
#include "Utils.hpp"
#endif

#ifndef INCLUDED_LSMPSa
#include "src/LSMPS/LSMPSa.hpp"
#endif

void saveResidual(Particle& par, int step){
    // Calculating the residual
    LSMPSa lsmpsa_du;
    LSMPSa lsmpsa_dv;
    lsmpsa_du.set_LSMPS(par.x, par.y, par.s, par.u, par.neighbor);
	lsmpsa_dv.set_LSMPS(par.x, par.y, par.s, par.v, par.neighbor);
    std::vector<double> _dudx = lsmpsa_du.get_ddx();
    std::vector<double> _dvdy = lsmpsa_dv.get_ddy();

    // Calculating vorticity
    std::vector<double> _dudy = lsmpsa_du.get_ddy();
    std::vector<double> _dvdx = lsmpsa_dv.get_ddx();

    std::vector<double> res(par.num,0.0);
    std::vector<double> vor(par.num,0.0);

    for (int i = 0; i < par.num; i++){
        res[i] = _dudx[i] + _dvdy[i];
        vor[i] = _dvdx[i] - _dudy[i];
    }

    // Saving the residual data
    std::ofstream data;
    std::string name,number;
    name  = "output/residual_";
    simUtil util_step;
    util_step.startCounter(step);
    util_step.saveName(number,step);
    name += number + ".csv";
    data.open(name);
    data << "x,y,vor,res\n";
    for (int i = 0; i < par.num; i++){
    data << ""  << par.x[i]
         << "," << par.y[i]
         << "," << vor[i]
         << "," << res[i]
         << "\n";
    }
    
    return;
}