#ifndef INCLUDED_INITIALIZATION
#define INCLUDED_INITIALIZATION
#include <fstream>
#include <string>
#include <sstream>

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

#ifndef INCLUDED_GEOMETRY
#include "../geometry/geometry.hpp"
#endif

// #ifndef INCLUDED_NEIGHBOR_BASE
// #include "../neighbor/base_grid.hpp"
// #endif

// #ifndef INCLUDED_NEIGHBOR
// #include "../neighbor/neighbor.hpp"
// #endif

class initialization
{
private:
    // Method to calculate particle interaction to geometry
    geometry d_geom;        // Package for geometry

    // 2D Particles initialization
    void init_2d_single_res(Particle &p);                               // [DONE]
    void init_2d_multi_res_single_block(const Body &b, Particle &p);    // [DONE]
    void init_2d_multi_res_multi_block(const Body &b, Particle &p);     // Ongoing
    void init_2d_multi_res_body_adjusted(const Body &b, Particle &p);   // [DONE]
    // void init_particle(const Body &b, Particle &p);                  // [TEST - Not Included]
    // void init_particle_test(Particle &p);                            // [TEST - Not Included]

    // 3D Particles initialization
    void init_domain_3d(Particle &p);

public:
    // Constructor De-constructor
    initialization(/* args */);
    ~initialization();
    
    // Initialization global procedure
    void generate(const Body &b, Particle &p);               // [DONE]
    void continue_simulation(const Body &b, Particle &p, int step);    // [NOT YET]
};
#endif