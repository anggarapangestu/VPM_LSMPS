#ifndef INCLUDED_INITIALIZATION
#define INCLUDED_INITIALIZATION

#include "../../Utils.hpp"
#include "../grid_block/gridNode.hpp"
#include "../geometry/geometry.hpp"
// #include "../remeshing/remeshing.hpp"

/**
 *  @brief  Initialization class consisted of method to generate 
 *  initial particle distribution.
 * 
 *  @headerfile initialization.hpp
*/
class initialization
{
private:
    // Private manager method
    void generate_particle(Particle &_particle, const std::vector<Body> &_bodyList, GridNode &_baseGrid);
    // void calculate_vorticity(Particle &_particle);

    // Particles initialization by grid block method
    void init_par_grid_block(Particle &_particle, const std::vector<Body> &_bodyList, GridNode &g);

    // **2D Particles initialization
    void init_2d_single_res(Particle &_particle);                                                   // [DONE]
    void init_2d_multi_res_single_block(Particle &_particle, const std::vector<Body> &_bodyList);    // [DONE]
    void init_2d_multi_res_multi_block(Particle &_particle, const std::vector<Body> &_bodyList);     // Ongoing
    void init_2d_multi_res_body_adjusted(Particle &_particle, const std::vector<Body> &_bodyList);   // [DONE]
    void init_2d_test_1(Particle &_particle);       // [Testing 1]
    void init_2d_test_2(Particle &_particle);       // [Testing 2]

    // **3D Particles initialization
    void init_3d_single_res(Particle &_particle);                                                    // [DONE]
    void init_3d_multi_res_single_block(Particle &_particle, const std::vector<Body> &_bodyList);    // [DONE]
    void init_3d_multi_res_body_adjusted(Particle &_particle, const std::vector<Body> &_bodyList);   // [DONE]

    // **Vorticity initialization type
    // 2D vorticity initial
    void perlman_vorticity(Particle &_particle);
    void eliptic_vorticity(Particle &_particle);
    void lamb_oseen_vorticity(Particle &_particle);
    // 3D vorticity initial
    void ring_vortex(Particle &_particle);
    void ring_vortex_oblique(Particle &_particle);

public:
    // **Initialization global procedure

    void calculate_vorticity(Particle &_particle);

    void initialize_particle(Particle &_particle, const std::vector<Body> &_bodyList, GridNode &_baseGrid);
    void initialize_vorticity(Particle &_particle, GridNode &_baseGrid, const std::vector<Body> &_bodyList);

    // **Update boundary particle
    void update_domain_boundary(Particle &_particle);

    // **Set active particle sign
    void set_active_sign(Particle &_particle, bool _flag);

    // **Particle reader
    void read_2d_particle(Particle &_particle, int _iteration);
    void read_3d_particle(Particle &_particle, int _iteration);

    // Debugging and testing subroutine
    // Solution to initial vorticity
    void perlman_velocity_solution(Particle &_particle, int _type);
    void lamb_oseen_solution(Particle &_particle, double _time, int _type);

    // A testing function
    void test_0(Particle &_particle);
    void test_1(Particle &_particle);
    void test_2(Particle &_particle);

    // A solution to laplace problem
    void initialize_laplace(Particle &_particle);
    void laplace(Particle &_particle);
};
#endif