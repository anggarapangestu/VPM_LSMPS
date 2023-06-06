# VPM_Boundary_Layer
Fluid particle base simulation for Boundary Layer on a flat plate.
The code is based on the previous work by Ical and Adhika.

## Program Scenario
Boundary layer is a steady flow phenomena as long the flat plate not moving (static).
This simulation using explicit approach in which simulating the boundary layer from it's formation.
Initially the domain is a fluid with uniform flow by property of ρ, μ, and Re.
The flat plate solid will be impulsively started in the domain then the boundary layer start to develop.
The wall boundary condition is simulated by Brinkman Penalization.

## Output
The output of this simulation is the steady state is 
- The Boundary Layer Formation
- Boundary Layer Velocity Profile
