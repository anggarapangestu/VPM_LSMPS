#ifndef INCLUDED_GEOMETRY
#define INCLUDED_GEOMETRY

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

class geometry
{
private:
	//2D Body Objects Geometry Generation
	void cylinder_generator(Body &b);		// [opt 1] Circular Cylinder
	void flat_plate_generator(Body &b);		// [opt 2] Flat Plate
	void normal_plate_generator(Body &b);	// [opt 3] Normal Plate
	void square_generator(Body &b);			// [opt 4] Square Cylinder
	void naca_generator(Body &b);			// [opt 5] NACA Airfoil
	
	//3D Body Objects Geometry Generation
	void ball_sphere(Body &b);				// [opt 1] Sphere

	// Other Package
	void u_var(Body &b);

public:
	// Geometry Method
	void generateBody(Body &b);
	void moving_body(int it, Body &b);

	// Minimum distance from particle to body surface
	void distance_calc(Particle &p, const Body& b);
	double distance_calc(std::vector<double> & pos, const Body& b);
};

#endif
