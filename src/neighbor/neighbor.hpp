#ifndef INCLUDED_NEIGHBOR
#define INCLUDED_NEIGHBOR

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

#ifndef INCLUDED_NEIGHBOR_BASE
#include "base_grid.hpp"
#endif

class neighbor
{
private:
	// The instances inside neighbor
	base_grid d_base_grid;
	CellList d_cell_list;

public:
	// The neighbor search manager
	void cell_list_init(Particle& parEval);
	void neigbor_search(Particle& parEval);
	bool particle_adaptation(const Particle& parEval, Particle& baseParticle, std::vector<double>& PARsize);
	std::vector<std::vector<int>> par2grid_neigbor_search(Particle& parEval, const Particle& parBase);

	//2-D Neighbour Search;
	// TODO: generate neighborhood using link-list algorithm
	std::vector<std::vector<int>> link_list(const int np, const std::vector<double> &sp,
											const std::vector<double> &xp, const std::vector<double> &yp, int neighbor_scale);
	// ! -- obsolete --
	//// void link_list(const int np, const std::vector<double> &sp, const std::vector<double> &xp, const std::vector<double> &yp,
	////     std::vector<int> &pair_i, std::vector<int> &pair_j, int neighbor_scale);
	// TODO: generate neighborhood using direct searching
	std::vector<std::vector<int>> direct_find(const int np, const std::vector<double> &sp,
											  const std::vector<double> &xp, const std::vector<double> &yp, const int neighbor_scale);
	// ! -- obsolete --
	//// void direct_find(const int np, const std::vector<double> &sp, const std::vector<double> &xp, const std::vector<double> &yp,
	//// 				 std::vector<int> &pair_i, std::vector<int> &pair_j, const int neighbor_scale);
	
	//3-D Neighbour Search
	//direct search
	void direct_3d(Particle &p, const int np, const double neighbor_scale);
	//Newest Method:Neighbour Search Using Spatial Hashing : by:Ical 
	void hash(Particle &p, Cell &c, const int np, double cell_size);
};

#endif
