#include "advection.hpp"

void advection::advection_euler(Particle &p)
{   
    // Advection prompt
    printf("\nCalculating Advection ...\n");

    // Advection computational time intialization
    clock_t t;
    t = clock();

    // Advection computation
    for (int i = 0; i < p.num; i++)
    {
        p.x[i] += p.u[i] * Pars::dt;
        p.y[i] += p.v[i] * Pars::dt;
    }
    
    // Display computational time
    t = clock() - t;
    printf("<-> Advection comp. time:              [%f s]\n", (double)t/CLOCKS_PER_SEC);
}


void advection::advection_rk2(Particle &p, std::vector<double> &dfdt)
{
    // Internal variables
    Particle _particle;
    std::vector<int> _index;

    printf("\nCalculating convection ...");
    clock_t t; // time variables
    t = clock();
    double u1, u2, u3, v1, v2, v3;
    for (int i = 0; i < p.num; i++)
    {
        u1 = (p.u[i] / 2);
        u2 = u1 / 2;
        u3 = u2 ;
        v1 = (p.v[i] / 2);
        v2 = v1 / 2;
        v3 = v3;

        p.x[i] += (p.u[i] + 2*u1 + 2*u2 + u3)*Pars::dt/6;
        p.y[i] += (p.v[i] + 2*v1 + 2*v2 + v3)*Pars::dt/6;
    }

    t = clock() - t;
    printf(" [%f s].\n", ((float)t) / CLOCKS_PER_SEC);
}