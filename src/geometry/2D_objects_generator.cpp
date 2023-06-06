#include "geometry.hpp"

/*/ ========================================== ///
Note in 2D geometry:
> Body node is arranged in positive convenction direction (CCW)
/// ========================================== /*/

// [opt 1] Circular Cylinder Geometry Generation
void geometry::cylinder_generator(Body &b)
{
	// internal variables
	double *_alfa = new double [Pars::n_a];
	double d_alfa = (360.0e0/Pars::n_a)*Pars::pi/180.0e0;
	double *min = new double[2];
	double *max = new double[2];

	// surface point list of position angle (theta)
	_alfa[0] = 0.0e0;
	for (int i = 0; i < Pars::n_a; i++)
	{
		_alfa[i+1] = _alfa[i] + d_alfa;
	}

	// resize x and y coordinate of obstacle
	b.x.clear(); b.x.resize(Pars::n_a);
	b.y.clear(); b.y.resize(Pars::n_a);
	b.x_n.clear(); b.x_n.resize(Pars::n_a);
	b.y_n.clear(); b.y_n.resize(Pars::n_a);
	b.x_m.clear(); b.x_m.resize(Pars::n_a);
	b.y_m.clear(); b.y_m.resize(Pars::n_a);
	b.num = Pars::n_a;
	b.n_panel = Pars::n_a;
	
	// Assign initial value
	min[0] = 0;  // minimum x position
	max[0] = 0;  // maximum x position
	min[1] = 0;  // minimum y position
	max[1] = 0;  // maximum y position

	// surface body point coordinate, panel midpoint coordinate, and panel normal vector
	for (size_t i = 0; i < Pars::n_a; i++)
	{
		b.x[i] = Pars::Df/2 * std::cos(_alfa[i]);
		b.y[i] = Pars::Df/2 * std::sin(_alfa[i]);

		// Evaluate extremes
		min[0] = min[0] < b.x[i] ? min[0] : b.x[i];  // minimum x position
		max[0] = max[0] > b.x[i] ? max[0] : b.x[i];  // maximum x position
		min[1] = min[1] < b.y[i] ? min[1] : b.y[i];  // minimum y position
		max[1] = max[1] > b.y[i] ? max[1] : b.y[i];  // maximum y position

		b.x_n[i] = std::cos(_alfa[i] + d_alfa/2);
		b.y_n[i] = std::sin(_alfa[i] + d_alfa/2);
		if (i > 0){
			b.x_m[i-1] = (b.x[i] + b.x[i-1])/2;
			b.y_m[i-1] = (b.y[i] + b.y[i-1])/2;
		}
	}
	b.x_m[Pars::n_a-1] = (b.x[Pars::n_a-1] + b.x[0])/2;
	b.y_m[Pars::n_a-1] = (b.y[Pars::n_a-1] + b.y[0])/2;
	
	// Assign the extremes point
	b.min_pos.emplace_back(min[0]);
	b.min_pos.emplace_back(min[1]);
	b.max_pos.emplace_back(max[0]);
	b.max_pos.emplace_back(max[1]);

	// deallocating memory
	delete[] _alfa, min, max;
}

// [opt 2] Flat Plate Geometry Generation
void geometry::flat_plate_generator(Body &b)
{
	// internal variables declaration
	double lenpanelx, lenpanely, L, H;
	double *x_max = new double[2];
    double *x_min = new double[2];
	int numpanelx, numpanely, N, id;

	// assigning internal variables
	L = Pars::Df;
	H = L*Pars::H_star;
	numpanelx = std::floor((L/(2*(L+H))) * Pars::n_a);
	numpanely = std::ceil(numpanelx*Pars::H_star);
	
	// calculating internal variables
	N = 2*(numpanelx + numpanely);
	lenpanelx = L/numpanelx;
	lenpanely = H/numpanely;
	// maximum vertex coordinate
		x_max[0] = L;			// x_max
		x_max[1] = H/2.0e0;		// y_max
	// minimum vertex coordinate
		x_min[0] = 0;			// x_min
		x_min[1] = -H/2.0e0;	// y_min

	// resize x and y coordinate of body object and body initialization
	b.x.clear(); b.x.resize(N);
	b.y.clear(); b.y.resize(N);
	b.x_m.clear(); b.x_m.resize(N);
	b.y_m.clear(); b.y_m.resize(N);
	b.x_n.clear(); b.x_n.resize(N);
	b.y_n.clear(); b.y_n.resize(N);
	b.num = N;
	b.n_panel = N;
	id = 0;		// The index of body node

	// surface body point coordinate, panel midpoint coordinate, and panel normal vector
	for (int i = 0; i < numpanelx; i++) // SEGMENT 1: Vertex node 1-2 (bottom segment)
	{
		b.x[id] = x_min[0] + i*lenpanelx;
		b.y[id] = x_min[1];
		b.x_m[id] = b.x[id] + lenpanelx/2;
		b.y_m[id] = x_min[1];
		b.x_n[id] = 0;
		b.y_n[id] = -1;
		id++;
	}
	for (int i = 0; i < numpanely; i++) // SEGMENT 2: Vertex node 2-3 (right segment)
	{
		b.x[id] = x_max[0];
		b.y[id] = x_min[1] + i*lenpanely;
		b.x_m[id] = x_max[0];
		b.y_m[id] = b.y[id] + lenpanely/2;
		b.x_n[id] = 1;
		b.y_n[id] = 0;
		id++;
	}
	for (int i = 0; i < numpanelx; i++) // SEGMENT 3: Vertex node 3-4 (top segment)
	{
		b.x[id] = x_max[0] - i*lenpanelx;
		b.y[id] = x_max[1];
		b.x_m[id] = b.x[id] - lenpanelx/2;
		b.y_m[id] = x_max[1];
		b.x_n[id] = 0;
		b.y_n[id] = 1;
		id++;
	}
	for (int i = 0; i < numpanely; i++) // SEGMENT 4: Vertex node 4-1 (left segment)
	{
		b.x[id] = x_min[0];
		b.y[id] = x_max[1] - i*lenpanely;
		b.x_m[id] = x_min[0];
		b.y_m[id] = b.y[id] - lenpanely/2;
		b.x_n[id] = -1;
		b.y_n[id] = 0;
		id++;
	}
	
	// Assign the extremes point
	b.min_pos.emplace_back(x_min[0]);
	b.min_pos.emplace_back(x_min[1]);
	b.max_pos.emplace_back(x_max[0]);
	b.max_pos.emplace_back(x_max[1]);

	// deallocating memory
	delete[] x_max;
	delete[] x_min;
}

// [opt 3] Normal to Freestream Plate Geometry Generation
void geometry::normal_plate_generator(Body &b)
{
	// internal variables declaration
	double lenpanelx, lenpanely, L, W;
	double *x_max = new double[2];
    double *x_min = new double[2];
	int numpanelx, numpanely, N, id;

	// assigning internal variables
	L = Pars::Df;
	W = L*Pars::H_star;
	numpanely = std::floor((L/(2*(L+W))) * Pars::n_a);
	numpanelx = std::ceil(numpanely*Pars::H_star);
	
	// calculating internal variables
	N = 2*(numpanelx + numpanely);
	lenpanely = L/numpanely;
	lenpanelx = W/numpanelx;
	// maximum vertex coordinate
		x_max[0] = W;			// x_max
		x_max[1] = L/2.0e0;		// y_max
	// minimum vertex coordinate
		x_min[0] = 0;			// x_min
		x_min[1] = -L/2.0e0;	// y_min

	// resize x and y coordinate of body object and body initialization
	b.x.clear(); b.x.resize(N);
	b.y.clear(); b.y.resize(N);
	b.x_m.clear(); b.x_m.resize(N);
	b.y_m.clear(); b.y_m.resize(N);
	b.x_n.clear(); b.x_n.resize(N);
	b.y_n.clear(); b.y_n.resize(N);
	b.num = N;
	b.n_panel = N;
	id = 0;		// The index of body node

	// surface body point coordinate, panel midpoint coordinate, and panel normal vector
	for (int i = 0; i < numpanelx; i++) // SEGMENT 1: Vertex node 1-2 (bottom segment)
	{
		b.x[id] = x_min[0] + i*lenpanelx;
		b.y[id] = x_min[1];
		b.x_m[id] = b.x[id] + lenpanelx/2;
		b.y_m[id] = x_min[1];
		b.x_n[id] = 0;
		b.y_n[id] = -1;
		id++;
	}
	for (int i = 0; i < numpanely; i++) // SEGMENT 2: Vertex node 2-3 (right segment)
	{
		b.x[id] = x_max[0];
		b.y[id] = x_min[1] + i*lenpanely;
		b.x_m[id] = x_max[0];
		b.y_m[id] = b.y[id] + lenpanely/2;
		b.x_n[id] = 1;
		b.y_n[id] = 0;
		id++;
	}
	for (int i = 0; i < numpanelx; i++) // SEGMENT 3: Vertex node 3-4 (top segment)
	{
		b.x[id] = x_max[0] - i*lenpanelx;
		b.y[id] = x_max[1];
		b.x_m[id] = b.x[id] - lenpanelx/2;
		b.y_m[id] = x_max[1];
		b.x_n[id] = 0;
		b.y_n[id] = 1;
		id++;
	}
	for (int i = 0; i < numpanely; i++) // SEGMENT 4: Vertex node 4-1 (left segment)
	{
		b.x[id] = x_min[0];
		b.y[id] = x_max[1] - i*lenpanely;
		b.x_m[id] = x_min[0];
		b.y_m[id] = b.y[id] - lenpanely/2;
		b.x_n[id] = -1;
		b.y_n[id] = 0;
		id++;
	}
	
	// Assign the extremes point
	b.min_pos.emplace_back(x_min[0]);
	b.min_pos.emplace_back(x_min[1]);
	b.max_pos.emplace_back(x_max[0]);
	b.max_pos.emplace_back(x_max[1]);

	// deallocating memory
	delete[] x_max;
	delete[] x_min;
}

// [opt 4] Square Cylinder Geometry Generation
void geometry::square_generator(Body &b)
{
	// internal variables declaration
	double lenpanel, L, H, pos_max, pos_min;
	int numpanel, N, id;

	// assigning internal variables
	L = Pars::Df;
	numpanel = Pars::n_a / 4;
	
	// calculating internal variables
	N = 4 * numpanel;
	lenpanel = L/numpanel;
	pos_max = L/2.0e0;	// maximum vertex coordinate
	pos_min = -L/2.0e0;	// minimum vertex coordinate

	// resize x and y coordinate of body object and body initialization
	b.x.clear(); b.x.resize(N);
	b.y.clear(); b.y.resize(N);
	b.x_m.clear(); b.x_m.resize(N);
	b.y_m.clear(); b.y_m.resize(N);
	b.x_n.clear(); b.x_n.resize(N);
	b.y_n.clear(); b.y_n.resize(N);
	b.num = N;
	b.n_panel = N;
	id = 0;		// The index of body node

	// surface body point coordinate, panel midpoint coordinate, and panel normal vector
	for (int i = 0; i < numpanel; i++) // SEGMENT 1: Vertex node 1-2 (bottom segment)
	{
		b.x[id] = pos_min + i*lenpanel;
		b.y[id] = pos_min;
		b.x_m[id] = pos_min + (i + 0.5e0) * lenpanel;
		b.y_m[id] = pos_min;
		b.x_n[id] = 0;
		b.y_n[id] = -1;
		id++;
	}
	for (int i = 0; i < numpanel; i++) // SEGMENT 2: Vertex node 2-3 (right segment)
	{
		b.x[id] = pos_max;
		b.y[id] = pos_min + i*lenpanel;
		b.x_m[id] = pos_max;
		b.y_m[id] = pos_min + (i + 0.5e0) * lenpanel;
		b.x_n[id] = 1;
		b.y_n[id] = 0;
		id++;
	}
	for (int i = 0; i < numpanel; i++) // SEGMENT 3: Vertex node 3-4 (top segment)
	{
		b.x[id] = pos_max - i*lenpanel;
		b.y[id] = pos_max;
		b.x_m[id] = pos_max - (i + 0.5e0) * lenpanel;
		b.y_m[id] = pos_max;
		b.x_n[id] = 0;
		b.y_n[id] = 1;
		id++;
	}
	for (int i = 0; i < numpanel; i++) // SEGMENT 4: Vertex node 4-1 (left segment)
	{
		b.x[id] = pos_min;
		b.y[id] = pos_max - i*lenpanel;
		b.x_m[id] = pos_min;
		b.y_m[id] = pos_max - (i + 0.5e0) * lenpanel;
		b.x_n[id] = -1;
		b.y_n[id] = 0;
		id++;
	}
	
	// Assign the extremes point
	b.min_pos.emplace_back(pos_min);
	b.min_pos.emplace_back(pos_min);
	b.max_pos.emplace_back(pos_max);
	b.max_pos.emplace_back(pos_max);
}

// [opt 5] NACA Airfoil Geometry Generation
void geometry::naca_generator(Body &b)
{
	double *x    = new double [Pars::n_a+1];
	double *y    = new double [Pars::n_a+1];
	double *zeta = new double [Pars::n_a+1];
	double *yt   = new double [Pars::n_a+1];
	double *yc   = new double [Pars::n_a+1];
	double *lc   = new double [Pars::n_a+1];
	double sin_theta, cos_theta, l_c, p_a1;
	double **xy_rot = new double *[Pars::n_a+1];
	double **xy     = new double *[Pars::n_a+1];
	
	for (size_t i = 0; i < Pars::n_a+1; i++)
	{
		xy_rot[i] = new double[2];
		xy[i] 		= new double[2];
	}
	double Rotz [2][2];	

	for (int i = 0; i < Pars::n_a+1; i++)
	{
		zeta[i] = 2.0e0 * Pars::pi / Pars::n_a * i;
		x[i]    = 0.5e0 * (std::cos(zeta[i]) + 1);
	}

	for (int i = Pars::n_a/2+1; i < Pars::n_a; i++)
	{
		if (x[i] < Pars::p_a)
		{
			if (Pars::p_a == 0.0e0)
			{
				yc[i] = 0.0e0;
			}
			else
			{
				p_a1  = Pars::p_a;
				yc[i] = Pars::m_a/std::pow(p_a1,2)*(2.0e0*p_a1*x[i] - std::pow(x[i],2));
			}
		}
		else
		{
			yc[i] = Pars::m_a/std::pow((1.0e0-Pars::p_a),2)*((1.0e0-2.0e0*Pars::p_a)+2.0e0*Pars::p_a*x[i]-std::pow(x[i],2));
		}

	    //yt[i] = Pars::t_a/0.2e0*(0.2969e0*sqrt(x[i]) - 0.126e0*x[i] - 0.35160e0*std::pow(x[i],2) + 0.2843e0*std::pow(x[i],3) - 0.1015e0*std::pow(x[i],4));
		yt[i] = Pars::t_a*5*(0.2969e0*sqrt(x[i]) - 0.126e0*x[i] - 0.35160e0*std::pow(x[i],2) + 0.2843e0*std::pow(x[i],3) - 0.1036e0*std::pow(x[i],4));
		l_c   = sqrt(std::pow((x[i]-x[i-1]),2) + std::pow((yc[i]-yc[i-1]),2));
	    sin_theta = (yc[i] - yc[i-1])/l_c;
	    cos_theta = (x[i] - x[i-1])/l_c;

	    x[i] = x[i] - yt[i]*sin_theta;
	    y[i] = yc[i] + yt[i]*cos_theta;
	    x[Pars::n_a-i] = x[i] + yt[i]*sin_theta;
	    y[Pars::n_a-i] = yc[i] - yt[i]*cos_theta;
	}

	y[0]     = 0;
	y[Pars::n_a]   = 0;
	y[Pars::n_a/2] = 0;
	x[0]     = 1;
	x[Pars::n_a]   = 1; // chord length
	x[Pars::n_a/2] = 0;

	// -- 2D Airfoil
	if(Pars::vib == 99){
	for (int i = 0; i < Pars::n_a+1; i++)
	{
		xy[i][0] = x[i];
		xy[i][1] = y[i];
	}
	}else{
	// -- 2D Airfoil, aerodynamic center pindah ke 0,0
	
	for (int i = 0; i < Pars::n_a+1; i++)
	{
		xy[i][0] = x[i] - Pars::acx;
		xy[i][1] = y[i] - Pars::acy;
	}

	}
	Rotz[0][0] =  std::cos(Pars::alpha_a*Pars::pi/180.0e0);
	Rotz[0][1] =  std::sin(Pars::alpha_a*Pars::pi/180.0e0);
	Rotz[1][0] = -std::sin(Pars::alpha_a*Pars::pi/180.0e0);
	Rotz[1][1] =  std::cos(Pars::alpha_a*Pars::pi/180.0e0);

	// resize x and y coordinate of obstacle
	b.x.clear(); b.x.resize(Pars::n_a + 1);
	b.y.clear(); b.y.resize(Pars::n_a + 1);
	
	for (int i = 0; i < Pars::n_a+1; i++)
	{
		xy_rot[i][0] = Rotz[0][0]*xy[i][0] + Rotz[0][1]*xy[i][1];
		xy_rot[i][1] = Rotz[1][0]*xy[i][0] + Rotz[1][1]*xy[i][1];
		
		b.x[i] = xy_rot[i][0];
		b.y[i] = xy_rot[i][1];
	}

	// -- clear internal variables
	delete[] x;
	delete[] y;
	delete[] zeta;
	delete[] yt;
	delete[] yc;
	delete[] lc;
	for (int i = 0; i < Pars::n_a+1; i++)
	{
		delete[] xy_rot[i];
		delete[] xy[i];
	}
	delete[] xy_rot, xy;
	// null
	x = NULL;
	y = NULL;
	zeta = NULL;
	yt  = NULL;
	yc  = NULL;
	lc  = NULL;
	xy_rot  = NULL; xy = NULL;

}

/*
Note in geometry creation
*/