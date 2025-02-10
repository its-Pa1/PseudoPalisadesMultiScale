
Main files/scripts:
The file "main_hyp_all.m" is the main matlab script which calls all the functions (scripts) and simulates the 2-D model on a square grid. It produces (most of) the results shown in the subsection 2.5.5 of "Diss_Kumar_Pawan.pdf" present in the parent directory.


Functions:
tissue_Q_macro.m: computes the macroscopic tissue Q described in subsection 2.5.2 of the thesis.
function_g.m: computes the g(S(x)) involved in eq. 2.4.46.
tissue_unimodal.m: computes the directed tissue (q) given by eq. 2.5.5
set_tumor_diff.m: assemble the tumor diffusion matrix and store other related variables
set_acidity_diff.m: assemble the tumor diffusion matrix.

Scripts:
main_SU_delta1.m: Simulates the model (2.4.46 & 2.3.7) with anisotropic tissue and with two occlusion sites as initial data and with delta = 1 in eq. 2.5.5, the results are shown in Figure 2.12 and 2.13.
main_SU_exp2_delta0p2.m: Simulates the model (2.4.46 & 2.3.7) with anisotropic tissue and with two occlusion sites as initial data and with delta = 0.2 in eq. 2.5.5, the results are shown in Figure 2.14 and 2.15.


Method: for time discretisation IMEX method is used. Diffusion is solved implicitly(implicit Euler) while source and taxis/advection terms explicitly(Euler).
 
For space discretisation, Weikert's 3X3 (central diff, 9 points as mentioned in the Table 2.3 of the thesis) stencil discretisation(for anisotropic diffusion) has been used for tumor diffusion and standard 5 point stencil(again central diff, as shown in Table 2.2)  for acid diffusion. Second order upwind scheme with Von Leer flux limiter is used in both x and y direction for the taxis term (see subsection 2.5.3 of the thesis).


Time: 30-35 minutes(on 8 core, 16gb RAM, i7-processor @ 3.60GHz speed, Linux) to 35-40 minutes(dual core, 8GB RAM, i5-processor, @2,3GHz speed, Mac)

Requirement: Matlab(2018a or later version, it should work for older version also)
