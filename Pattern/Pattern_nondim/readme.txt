Main files/scripts:
The file "Pattern_main_i.m" (i=const, D1, D2)is the main matlab script which calls all the functions and simulates the non dim 1-D model for different diffusion coefficients D (const, D1 and D2 as mentioned in subsection 2.5.6 of the thesis) and plot tumor & acidity w.t.to time and space

Functions:
set_diff_mat_const.m: assemble the tumor diffusion matrix and store const D
set_diff_mat_D1.m: assemble the tumor diffusion matrix and store D1
set_diff_mat_D2.m: assemble the tumor diffusion matrix and store D2

Method: for time discretisation IMEX method is used. Diffusion is solved implicitly(implicit Euler) while source and taxis/advection terms explicitly(Euler).
 
For space discretisation, standard 5 point stencil(central diff)  for diffusion is used
Upwind scheme(first order) is used for the taxis term.


Requirement: Matlab(2018a or later version, it should work for older version also)