Main files/scripts:
The file "main_diff.m" is the main matlab script which calls "main_2D.m" and "main_2D_without_g.m", these two scripts simulate the model (2.4.16 & 2.3.7) with and without g(S) respectively. The produced results are shown in Figure 2.11 of the thesis.



Functions:
tissue_Q_macro.m: computes the macroscopic tissue Q described in subsection 2.5.2 of the thesis.
function_g.m: computes the g(S(x)) involved in eq. 2.4.16.
tissue_q_un.m: computes the undirected tissue small q gievn by eq. 2.5.4
set_tumor_diff.m: assemble the tumor diffusion matrix and store other related variables
set_acidity_diff.m: assemble the tumor diffusion matrix.

Method: for time discretisation IMEX method is used. Diffusion is solved implicitly(implicit Euler) while source and taxis/advection terms explicitly(Euler).
 
For space discretisation, Weikert's 3X3 (central diff, 9 points as mentioned in the Table 2.3 of the thesis) stencil discretisation(for anisotropic diffusion) has been used for tumor diffusion and standard 5 point stencil(again central diff, as shown in Table 2.2)  for acid diffusion. Upwind scheme(first order) is used in both x and y direction for the taxis term (see subsection 2.5.3 of the thesis).

Time: 8-10 minutes(on 8 core, 16gb RAM, i7-processor @ 3.60GHz speed, Linux) to 10-12 minutes(dual core, 8GB RAM, i5-processor, @2,3GHz speed, Mac)

Requirement: Matlab(2018a or later version, it should work for older version also)