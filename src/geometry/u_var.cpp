#include "geometry.hpp"

/* ===== Code Description ===== //
This code assign the body velocity for each time step.
Body velocity values is given in the global.cpp parameter 
that are inputed manually in the code
*/

void geometry::u_var(Body &b)
{
	// Internal variables
	double *t = new double[Pars::nt];
	double t1, t2, TT, V_inf;

	printf("Calculating body velocity ... \n");
	clock_t d_time = clock();

	// Set variable's size
	b.uT.clear();
	b.uT.resize(Pars::nt);
	b.vT.clear();
	b.vT.resize(Pars::nt);

	// Body velocity in x-direction
	if (Pars::opt_u == 1)		// First type: Sudden move then stop
	{
		for (size_t i = 0; i < Pars::nt; i++)
		{
			t[i] = Pars::dt * (i + 1);
			if (t[i] <= 1.0e0)
				b.uT[i] = Pars::ubody;
			else // the body motion stop at t=1s until simulation finish
				b.uT[i] = 0.0e0;
		}
	}
	else if (Pars::opt_u == 2)	// Second type: Smooth motion
	{
		V_inf = Pars::ubody;
		t1 = 0.1e0;
		t2 = 1.1e0;
		TT = 0.01e0;

		for (size_t i = 0; i < Pars::nt; i++)
		{
			t[i] = Pars::dt * (i + 1);

			if ((t[i] <= (t1 + t2) * 0.5e0) && ((t[i] > t1)))
				b.uT[i] = V_inf * (1.0e0 - std::exp(-(t[i] - t1) / TT));
			else if ((t[i] >= (t1 + t2) * 0.5e0) && ((t[i] < t2)))
				b.uT[i] = V_inf * (1.0e0 - std::exp((t[i] - t2) / TT));

			if ((t[i] <= t1))
				b.uT[i] = V_inf * (std::exp((t[i] - t1) / TT)) / 2.0e0;
			else if ((t[i] >= t2))
				b.uT[i] = V_inf * (std::exp(-(t[i] - t2) / TT)) / 2.0e0;
		}
	}
	else if (Pars::opt_u == 3)	// Third type: Constant velocity
	{
		for (size_t i = 0; i < Pars::nt; i++)
			b.uT[i] = Pars::ubody;
	}
	
	// Body velocity in y-direction
	for (size_t i = 0; i < Pars::nt; i++)
	{
		b.vT[i] = Pars::vbody;	// Only constant velocity
	}

	d_time = clock() - d_time; // [ND]
	printf("<-> body velocity calculation\n");
	printf("    comp. time:                        [%f s]\n", (double)d_time / CLOCKS_PER_SEC);
	
	delete[] t;
}
