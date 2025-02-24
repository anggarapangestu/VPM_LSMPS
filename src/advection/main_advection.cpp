#include "advection.hpp"
#include "../LSMPS/LSMPSb.hpp"
#include "../neighbor/neighbor.hpp"

#define ADVECTION_TYPE 1    // The type of advection calculation
#define PC_or_RK 0     // Select the type for 0:=Predictor-Corrector or 1:=Runge-Kutta

// The basic advection
void advection::main_advection(Particle &p){
    // Advection prompt
    printf("\nCalculating advection ...\n");

    // Advection computational time manager
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    if (ADVECTION_TYPE == 1){
        printf("%s<+> Type 1: Euler type %s\n", FONT_CYAN, FONT_RESET);
        this->advection_euler(p);
    }
    else if (ADVECTION_TYPE == 2){
        printf("%s<+> Type 2: Runge Kutta 2nd order type %s\n", FONT_CYAN, FONT_RESET);
        this->advection_rk2(p);
    }
    
    // Display computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Advection computational time:      [%f s]\n", _time);

    return;
}

void advection::main_advection_corr(Particle &pPred, Particle &p){
    // Advection prompt
    printf("\nCalculating advection ...\n");

    // Advection computational time manager
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // MAIN CODE:
    // **********
    // This is the main code for 
    // Advection computation 2D
    if (DIM == 2){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        #if (PC_or_RK == 0)
            // Original
            p.x[i] = pPred.x[i] + 0.5 * (pPred.u[i] + p.u[i]) * Pars::dt;
            p.y[i] = pPred.y[i] + 0.5 * (pPred.v[i] + p.v[i]) * Pars::dt;
        #elif (PC_or_RK == 1)
            // For half step RK Test
            p.x[i] = pPred.x[i] + p.u[i] * Pars::dt;
            p.y[i] = pPred.y[i] + p.v[i] * Pars::dt;
        #endif
    }}
    // Advection computation 3D
    else if (DIM == 3){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        #if (PC_or_RK == 0)
            // Original
            p.x[i] = pPred.x[i] + 0.5 * (pPred.u[i] + p.u[i]) * Pars::dt;
            p.y[i] = pPred.y[i] + 0.5 * (pPred.v[i] + p.v[i]) * Pars::dt;
            p.z[i] = pPred.z[i] + 0.5 * (pPred.w[i] + p.w[i]) * Pars::dt;
        #elif (PC_or_RK == 1)
            // For half step RK Test
            p.x[i] = pPred.x[i] + p.u[i] * Pars::dt;
            p.y[i] = pPred.y[i] + p.v[i] * Pars::dt;
            p.z[i] = pPred.z[i] + p.w[i] * Pars::dt;
        #endif
    }}

    // Display computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Advection computational time:      [%f s]\n", _time);
    return;
}

void advection::main_advection_AB(Particle &pPrev, Particle &p){
    // Advection prompt
    printf("\nCalculating advection ...\n");

    // Advection computational time manager
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // MAIN CODE:
    // **********
    // This is the main code for 
    // Advection computation 2D
    if (DIM == 2){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        p.x[i] += (1.5*p.u[i] - 0.5*pPrev.u[i]) * Pars::dt;
        p.y[i] += (1.5*p.v[i] - 0.5*pPrev.v[i]) * Pars::dt;
    }}
    // Advection computation 3D
    else if (DIM == 3){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        p.x[i] += (1.5*p.u[i] - 0.5*pPrev.u[i]) * Pars::dt;
        p.y[i] += (1.5*p.v[i] - 0.5*pPrev.v[i]) * Pars::dt;
        p.z[i] += (1.5*p.w[i] - 0.5*pPrev.w[i]) * Pars::dt;
    }}

    // Display computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Advection computational time:      [%f s]\n", _time);
    return;
}

// The euler advection
// The velocity is constant and the particle is moved
void advection::advection_euler(Particle &p)
{   
    // Advection computation 2D
    if (DIM == 2){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        #if (PC_or_RK == 0)
            // Original PC
            p.x[i] += p.u[i] * Pars::dt;
            p.y[i] += p.v[i] * Pars::dt;
        #elif (PC_or_RK == 1)
            // Only for RK half step
            p.x[i] += p.u[i] * Pars::dt/2.0;
            p.y[i] += p.v[i] * Pars::dt/2.0;
        #endif
    }}
    // Advection computation 3D
    else if (DIM == 3){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        #if (PC_or_RK == 0)
            // Original PC
            p.x[i] += p.u[i] * Pars::dt;
            p.y[i] += p.v[i] * Pars::dt;
            p.z[i] += p.w[i] * Pars::dt;
        #elif (PC_or_RK == 1)
            // Only for RK half step
            p.x[i] += p.u[i] * Pars::dt/2.0;
            p.y[i] += p.v[i] * Pars::dt/2.0;
            p.z[i] += p.w[i] * Pars::dt/2.0;
        #endif
    }}
    
    return;
}

// Advection runge kutta 2nd order method
void advection::advection_rk2(Particle &p)
{   
    // Advection computation 2D
    if (DIM == 2){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        double u1 = (p.u[i] / 2);
        double u2 = u1 / 2;
        double u3 = u2 ;
        double v1 = (p.v[i] / 2);
        double v2 = v1 / 2;
        double v3 = v2;

        p.x[i] += (p.u[i] + 2*u1 + 2*u2 + u3)*Pars::dt/6;
        p.y[i] += (p.v[i] + 2*v1 + 2*v2 + v3)*Pars::dt/6;
    }}
    
    // Advection computation 3D
    else if (DIM == 3){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        double u1 = (p.u[i] / 2);
        double u2 = u1 / 2;
        double u3 = u2 ;
        double v1 = (p.v[i] / 2);
        double v2 = v1 / 2;
        double v3 = v2;
        double w1 = (p.w[i] / 2);
        double w2 = w1 / 2;
        double w3 = w2;

        p.x[i] += (p.u[i] + 2*u1 + 2*u2 + u3)*Pars::dt/6;
        p.y[i] += (p.v[i] + 2*v1 + 2*v2 + v3)*Pars::dt/6;
        p.z[i] += (p.w[i] + 2*w1 + 2*w2 + w3)*Pars::dt/6;
    }}
    return;
}

/**
 * @brief This is a particle point velocity data interpolation.
 *        It works by interpolating the velocity from the previous 
 *        field before advection, into the particle distribution after 
 *        advection.
 * @param _particleInit  The particle distribution before advection
 * @param _particleAdv   The particle distribution after advection [target interpolation]
 */
void advection::intplt_velocity(Particle &_particleInit, Particle &_particleAdv, const GridNode &_baseGrid){
    // Velocity Interpolation prompt
    printf("\nInterpolating velocities ...\n");
    
    // Copy the code from redistribution
    // Set up the neighbor for data interpolation
    std::vector<std::vector<int>> sourceNghList; // The neighbor of base particle using the ID of new particle (inter neighbor)
    std::vector<std::vector<int>> selfNghList;   // The neighbor of base particle relative to itself (self neihgbor)

    // Create an alias for each SOURCE data
    // Coordinates
    std::vector<double> &_srcX = _particleInit.x;
    std::vector<double> &_srcY = _particleInit.y;
    std::vector<double> &_srcZ = _particleInit.z;
    std::vector<double> &_srcSize = _particleInit.s;
    // Properties
    std::vector<double> &_srcU = _particleInit.u;
    std::vector<double> &_srcV = _particleInit.v;
    std::vector<double> &_srcW = _particleInit.w;   // 3D Only

    // Create an alias for each TARGET data
    // Coordinates
    std::vector<double> &_tarX = _particleAdv.x;
    std::vector<double> &_tarY = _particleAdv.y;
    std::vector<double> &_tarZ = _particleAdv.z;
    std::vector<double> _tarSize;                   // = _particleAdv.s;
    // Properties
    std::vector<double> &_tarU = _particleAdv.u;
    std::vector<double> &_tarV = _particleAdv.v;
    std::vector<double> &_tarW = _particleAdv.w;    // 3D Only

    // NOTE:
    //  > Both [self neighbor] and [inter neighbor] hold the same value if it is not adapted
    //  > If adaptation is performed, the "Base Particle" hold the [inter neighbor] data, while 
    //     [self neighbor] need to be calculated next
    //  > The [inter neighbor] is different if adapted and not adapted
    //     If adapted it already obtain all

    // Here pretend that the neighbor is just the same (Must be maintained by courant less than 0.5)
    // sourceNghList = _particleAdv.neighbor;
    selfNghList = _particleAdv.neighbor;
    MESSAGE_LOG << "Evaluating Inter Neighbor\n";
    if (Pars::opt_neighbor == 4){
        neighbor nghTools;
        nghTools.inter_neigbor_search(_particleInit, _particleAdv, _baseGrid, 1, sourceNghList, _tarSize);
        // _particleAdv.s = _tarSize;
        // sourceNghList = _particleAdv.neighbor;
    }else{
        neighbor nghTools;
        nghTools.neigbor_search(_particleAdv, _particleInit, sourceNghList);
        _tarSize = _particleAdv.s;
    }

    // Add the particle ID itself into the neighborID
    if (!Pars::flag_ngh_include_self){
        // ERROR_LOG << "GETTING INTO THIRD PAGE\n";
        for (int i = 0; i < _particleAdv.num; i++){
            selfNghList[i].push_back(i);
            // sourceNghList[i].push_back(i);
        }
    }

    // std::ofstream _write;
    // _write.open("ngh.csv");
    // int __ID = 23775;
    // _write << "x,y,s\n";
    // _write << _particleAdv.x.at(__ID) << ","
    //        << _particleAdv.y.at(__ID) << ","
    //        << _particleAdv.s.at(__ID) << ","
    //        << "\n";
    // for (int i = 0; i < sourceNghList.at(__ID).size(); i++){
    // int& _nghID = sourceNghList.at(__ID).at(i);
    // _write << _particleInit.x.at(_nghID) << ","
    //        << _particleInit.y.at(_nghID) << ","
    //        << _particleInit.s.at(_nghID) << ","
    //        << "\n";
    // }
    // _write.close();

    // std::ofstream _write;
    // _write.open("ngh_.csv");
    // __ID = 2255;
    // _write << "x,y,s\n";
    // _write << _particleAdv.x.at(__ID) << ","
    //        << _particleAdv.y.at(__ID) << ","
    //        << _particleAdv.s.at(__ID) << ","
    //        << "\n";
    // for (int i = 0; i < sourceNghList.at(__ID).size(); i++){
    // int& _nghID = sourceNghList.at(__ID).at(i);
    // _write << _particleInit.x.at(_nghID) << ","
    //        << _particleInit.y.at(_nghID) << ","
    //        << _particleInit.s.at(_nghID) << ","
    //        << "\n";
    // }
    // _write.close();
    // exit(0);

    // std::ofstream _write;
    // _write.open("ngh.csv");
    // int __ID = 1507;        // The id of grid particle
    // _write << "x,y,s\n";
    // _write << _particleInit.x.at(__ID) << ","
    //        << _particleInit.y.at(__ID) << ","
    //        << _particleInit.s.at(__ID) << ","
    //        << "\n";
    // for (int i = 0; i < sourceNghList.at(__ID).size(); i++){
    // int& _nghID = sourceNghList.at(__ID).at(i);
    // _write << _particleAdv.x.at(_nghID) << ","
    //        << _particleAdv.y.at(_nghID) << ","
    //        << _particleAdv.s.at(_nghID) << ","
    //        << "\n";
    // }
    // _write.close();

    // _write.open("ngh_.csv");
    // __ID = 34428;        // The id of grid particle
    // _write << "x,y,s\n";
    // _write << _particleInit.x.at(__ID) << ","
    //        << _particleInit.y.at(__ID) << ","
    //        << _particleInit.s.at(__ID) << ","
    //        << "\n";
    // for (int i = 0; i < sourceNghList.at(__ID).size(); i++){
    // int& _nghID = sourceNghList.at(__ID).at(i);
    // _write << _particleAdv.x.at(_nghID) << ","
    //        << _particleAdv.y.at(_nghID) << ","
    //        << _particleAdv.s.at(_nghID) << ","
    //        << "\n";
    // }
    // _write.close();
    // exit(0);

    // PROCEDURE 2: Interpolation of Velocity
    // ************
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // Perform LSMPS interpolation of vorticity
    LSMPSb lsmpsb;
    if (DIM == 2){
        // Interpolate the x velocity
        lsmpsb.set_LSMPS(_tarX, _tarY, _tarSize, _tarU,
                         _srcX, _srcY, _srcSize, _srcU,
                         sourceNghList, selfNghList);
        // Update the particle base data vorticity and vortex strength
        _particleAdv.u = lsmpsb.get_d0();

        // Interpolate the y velocity
        lsmpsb.set_LSMPS(_tarX, _tarY, _tarSize, _tarV,
                         _srcX, _srcY, _srcSize, _srcV,
                         sourceNghList, selfNghList);
        // Update the particle base data vorticity and vortex strength
        _particleAdv.v = lsmpsb.get_d0();
    }
    else if (DIM == 3){
        // Interpolate the velocity in x direction
        lsmpsb.set_LSMPS_3D(_tarX, _tarY, _tarZ, _tarSize, _tarU,
                            _srcX, _srcY, _srcZ, _srcSize, _srcU,
                            sourceNghList, selfNghList);
        _particleAdv.u = lsmpsb.get_d0();

        // Interpolate the velocity in y direction
        lsmpsb.set_LSMPS_3D(_tarX, _tarY, _tarZ, _tarSize, _tarV,
                            _srcX, _srcY, _srcZ, _srcSize, _srcV,
                            sourceNghList, selfNghList);
        _particleAdv.v = lsmpsb.get_d0();

        // Interpolate the velocity in y direction
        lsmpsb.set_LSMPS_3D(_tarX, _tarY, _tarZ, _tarSize, _tarW,
                            _srcX, _srcY, _srcZ, _srcSize, _srcW,
                            sourceNghList, selfNghList);
        _particleAdv.w = lsmpsb.get_d0();
    }
    
    // Time manager console log
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Calculating LSMPS interpolation:   [%f s]\n", _time);

    return;
}