#include "geometry.hpp"

/**
 *  @brief Generate all body data based on body object parameter in "global.hpp".
 *  NOTE: Not affected by N_BODY.
 *  
 *  @param	_bodyList	The list of body data container to be generated.
 * 
 *  @headerfile geometry.hpp
 */
void geometry::generateBody(std::vector<Body> &_bodyList){
    printf("+------------------ Body Region ------------------+\n");	// Opening section

    // Create all body part
	for (size_t i = 0; i < _bodyList.size(); i++){
		// Aliasing variable
		int part = i;
		Body &b = _bodyList.at(part);
		
		MESSAGE_LOG;
		printf("Generating surface node of body %d ... \n", part);
		
		// Computational time accumulation
		#if (TIMER_PAR == 0)
			// Timer using super clock
			std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
			// // Timer using simple timer
			// clock_t tick = clock();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			double _time = omp_get_wtime();
		#endif

		
		// Create the solid body data according to body option in Pars namespace
		if (DIM == 2)
		{
			// Generate 2 Dimensional Body Geometry
			switch (Pars::opt_body[part]){
			case 1:	// Circular cylinder
				printf("%sType 1: Circular cylinder %s\n", FONT_CYAN, FONT_RESET);
				this->cylinder_generator(b, part);
				break;
			case 2:	// Square cylinder
				printf("%sType 2: Square cylinder %s\n", FONT_CYAN, FONT_RESET);
				this->square_generator(b, part);
				break;
			case 3:	// Normal plate
				printf("%sType 3: Normal plate %s\n", FONT_CYAN, FONT_RESET);
				this->normal_plate_generator(b, part);
				break;
			case 4:	// Flat plate
				printf("%sType 4: Flat plate %s\n", FONT_CYAN, FONT_RESET);
				this->flat_plate_generator(b, part);
				break;
			case 5:	// NACA 4-series airfoil
				printf("%sType 5: NACA-%d%d%d airfoil %s\n", FONT_CYAN, 
									int(100*Pars::m_a),
									int( 10*Pars::p_a),
									int(100*Pars::t_a), FONT_RESET);
				this->naca_generator(b, part);
				break;
			default:
				break;
			}			
		}
		else if (DIM == 3)
		{
			// Generate 3 Dimensional Body Geometry
			switch (Pars::opt_body[part]){
			case 1:	// Sphere
				printf("%sType 1: Sphere %s\n", FONT_CYAN, FONT_RESET);
				this->ball_sphere_generator(b, part);
				break;
			case 2:	// Cube
				printf("%sType 2: Cube %s\n", FONT_CYAN, FONT_RESET);
				this->cube_generator(b, part);
				break;
			case 3: // 3D normal plate
				printf("%sType 3: 3D normal plate %s\n", FONT_CYAN, FONT_RESET);
				this->normal_plate_3d_generator(b, part);	// To be create further
				break;
			case 4: // 3D flat plate
				printf("%sType 4: 3D flat plate %s\n", FONT_CYAN, FONT_RESET);
				this->flat_plate_3d_generator(b, part);		// To be create further
				break;
			case 5:	// Torus
				printf("%sType 5: Torus %s\n", FONT_CYAN, FONT_RESET);
				this->torus_generator(b, part);
				break;
			case 6: // 3D Heart
				printf("%sType 6: Heart %s\n", FONT_CYAN, FONT_RESET);
				this->heart_generator(b, part);
				break;
			default:
				break;
			}
		}
		
		// Geometry computational time display
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
			double _time = span.count();
			// // Timer using simple clock
			// double _time = (clock() - tick)/CLOCKS_PER_SEC;
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime() - _time;
		#endif
		printf("<-> Solid body generation comp. time:  [%f s]\n", _time);

		// Update the geometry velocity [OLD package but correlated in penalization]
		this->define_vel_variation(b);
		printf("\n");
	}

    printf("+-------------------------------------------------+\n\n");	// Closing section
}

/*
Note for geometry
*/
