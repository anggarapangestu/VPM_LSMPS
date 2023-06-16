#ifndef INCLUDED_DATA_SAVING
#define INCLUDED_DATA_SAVING
#include <iomanip> // std::setw, std::setfill
#include <sstream> // std::stringstream
#include <fstream> // std::c_str()

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

#ifndef INCLUDED_LSMPS
#include "../LSMPS/LSMPSa.hpp"
#endif

#ifndef INCLUDED_BIOT_SAVART
#include "../velocity_poisson/velocity_biot_savart.hpp"
#endif

#ifndef INCLUDED_NEIGHBOR
#include "../neighbor/neighbor.hpp"
#endif

#ifndef INCLUDED_PENALIZATION
#include "../penalization/penalization.hpp"
#endif

#ifndef INCLUDED_SAVE_DATA_BASE
#include "base_save_data.hpp"
#endif

using namespace std;

class neighbor;
class base_save_data;
class save_data
{
#define w10 std::setw(10) // spare width while saving data
#define w16 std::setw(16) // spare width while saving data
#define w20 std::setw(20) // spare width while saving data

#pragma region instances
	velocity_biot_savart d_base_poisson;
	base_grid d_base_grid;
	// base_remeshing d_base_remeshing;
	penalization d_penalization;
	neighbor d_neighbor;
	base_save_data d_base_save_data;
#pragma endregion

#pragma region internal_variables
	// internal variables
	// @param; lifetime = singleton
	

	// ...variables...

	// @param; lifetime = transient
	
	std::vector<double> _Isumx;
	std::vector<double> _Isumy;

	std::vector<int> a1a2;
	std::vector<int> a2a3;
	std::vector<int> a3a4;
	std::vector<int> a4a1;
	std::vector<int> particles_inside;

	std::vector<double> u_jmin1;                                
	std::vector<double> u_jmax1;	
	std::vector<double> v_imin1;                                
	std::vector<double> v_imax1;	

	std::vector<double> u_jmin2;                                
	std::vector<double> u_jmax2;	
	std::vector<double> v_imin2;                                
	std::vector<double> v_imax2;

	std::vector<double> _omegasumx;
	std::vector<double> _omegasumy;

	std::vector<int> point_lists;

	Particle _p;

#pragma endregion

	// -- local functions

	void force_linear_impulse(int it, int np, const vector<double> &xp, const vector<double> &yp,
							  const vector<double> &gpz, const vector<double> &sp);

	// void grid_data(int np, vector<double> &xp, vector<double> &yp, vector<double> &sp, vector<double> &gpz,
	// 			   vector<double> &up, vector<double> &vp, int &nvertn, vector<vector<double>> &vertex);

public:
	// Basic saving data
	void save_state(Particle &p, string s, int type);		// Saving the state of particle distribution
	void save_state(const Body &b, string s);				// Saving the state of geometry
	void summary_log();
	void save_summary_log(Particle& par);
	
	//void particle_data_reading(int it, int &np, vector<double> &xp, vector<double> &yp, vector<double> &sp, vector<double> &gpz, vector<double> &up, vector<double> &vp);
	//void particle_data_reading(int &nvertn, std::vector<std::vector<double>> &vertex);
	
	void output(int it, Particle &p, Body &b, double &cum_time);

	void Force2(int iT, double A1, double A2, double A3, double A4, int nx, int ny, Particle &p);
	void set_variables(int numberofparticles_x, int numberofparticles_y, int totalnumberofparticles);
	// void save_dat(Particle &p, int np, vector<int> &numbar, vector<int> &batas_bawah, vector<int> &batas_atas, vector<int> &batas_kiri, vector<int> &batas_kanan, string s);
	void insert_data(Particle &p, int i);
	void sep(Particle p, vector<double> a, vector<double> b, vector<double> c, vector<double> d, vector<double> e, string s);
	void save_dat(Particle &p, int np, vector<int> &numbar, vector<int> &batas_bawah, vector<int> &batas_atas, 
				  vector<int> &batas_kiri, vector<int> &batas_kanan, string s, vector<double> &satu, vector<double> &dua,
				  vector<double> &tiga, vector<double> &empat, vector<double> &lima );

};

#endif
