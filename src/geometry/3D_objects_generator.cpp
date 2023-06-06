#include "geometry.hpp"

void geometry::ball_sphere(Body &b){
    // internal variables
	double *_x    = new double [Pars::n_a+1];
	double *_y    = new double [Pars::n_a+1];
    double *_z    = new double [Pars::n_a+1];
	double *_alfa = new double [Pars::n_a+1];

	// resize x and y coordinate of obstacle
	b.x.clear(); //b.x.resize(Pars::n_a + 1);
	b.y.clear(); //b.y.resize(Pars::n_a + 1);

	//#pragma omp parallel for
	/*for (size_t i = 0; i < Pars::n_a+1; i++)
	{
		// cylinder with a unit diameter
		_x[i] = Pars::lx/2 * std::cos(_alfa[i]);
		_y[i] = Pars::ly/2 * std::sin(_alfa[i]);
		b.x[i] = _x[i];
		b.y[i] = _y[i];
	}*/
    int bruh = 0;
    for(double phi = 0.; phi < 2* Pars::pi; phi += Pars::pi/30.){ //azimuth [0, 2pi]
        for(double theta = 0.; theta < Pars::pi; theta += Pars::pi/30.)//elevation [0, pi]
        {
            double x =  Pars::Df/2 * std::cos(phi) * std::sin(theta);
            double y =  Pars::Df/2 * std::sin(phi) * std::sin(theta);
            double z =  Pars::Df/2 * std::cos(theta);
            b.x.push_back(x);
            b.y.push_back(y);
            b.z.push_back(z);
            bruh++;
        } 
    }

    // deallocating memory
	delete[] _x;
	delete[] _y;
	delete[] _alfa;

	// // saving data
    // std::string nama;
	// std::ofstream ofs;

	// nama.append("3D_body");
	// nama.append(".csv");
    // ofs.open(nama.c_str());
    // ofs << "" << "i" << ","  << "xp" << "," << "yp" <<","<< "zp\n" ;
    // for (int i = 0; i < bruh; i++)
    // {
    //     ofs << "" << i
    //         << "," << b.x[i]
    //         << "," << b.y[i]
    //         << "," << b.z[i]<<"\n";
    // }
        
    // ofs.close();
}

/*
//counter clockwise (?)
void make_Sphere(VERTEX center, double r, std::vector<VERTEX> &spherePoints)
{
    const double PI = 3.141592653589793238462643383279502884197;
    spherePoints.clear();

    // Iterate through phi, theta then convert r,theta,phi to  XYZ
    for (double phi = 0.; phi < 2*PI; phi += PI/10.) // Azimuth [0, 2PI]
    {
        for (double theta = 0.; theta < PI; theta += PI/10.) // Elevation [0, PI]
        {
            VERTEX point;
            point.x = r * cos(phi) * sin(theta) + center.x;
            point.y = r * sin(phi) * sin(theta) + center.y;
            point.z = r            * cos(theta) + center.z;
            spherePoints.push_back(point);        
        }
    }
    return;
}

void geometry::pipe(Body &b){
}
*/