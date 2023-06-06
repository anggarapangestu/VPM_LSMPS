#include "velocity_biot_savart.hpp"
#include <omp.h>

void velocity_biot_savart::biotsavart_direct_2d(Particle &pi, Particle &pj)
{
	// internal variables
	double dx, dy;
	double rij2, rij, sij, sij2;
	double gsig;

	ProgressBar progressBar(pi.num, 35, '|', ' ');

	//#pragma omp parallel for
	for (int i = 1; i <= pi.num; i++)
	
	{
		++progressBar;		   // record the tick
		progressBar.display(); // display the bar

		for (int j = 1; j <= pj.num; j++)
		{
			dx = pi.x[i - 1] - pj.x[j - 1];
			dy = pi.y[i - 1] - pj.y[j - 1];
			rij2 = std::pow(dx, 2) + std::pow(dy, 2);
			rij = std::sqrt(rij2);

			if (rij2 > 0.0e0)
			{
				sij2 = (std::pow(pi.s[i - 1], 2) + std::pow(pj.s[j - 1], 2)) * 0.25e0;
				sij = std::sqrt(sij2);
				regul_func_2d(rij, sij, gsig);

				pi.u[i - 1] = pi.u[i - 1] - (dy * pj.gz[j - 1]) * gsig / rij2;
				pi.v[i - 1] = pi.v[i - 1] - (-dx * pj.gz[j - 1]) * gsig / rij2;
			}
		}
	}
	progressBar.done();
}
