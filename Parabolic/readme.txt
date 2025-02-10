Main files/scripts:
The file "main_all.m" is the main matlab script which calls all the functions and simulates the 2-D model on a square grid. It produces (most of) the results shown in the subsection 2.5.4 of "Diss_Kumar_Pawan.pdf" present in the parent directory.


Functions:
tissue_Q_macro.m: computes the macroscopic tissue Q described in subsection 2.5.2 of the thesis.
function_g.m: computes the g(S(x)) involved in eq. 2.4.16.
tissue_q_un.m: computes the undirected tissue small q gievn by eq. 2.5.4
set_tumor_diff.m: assemble the tumor diffusion matrix and store other related variables
set_acidity_diff.m: assemble the tumor diffusion matrix

Scripts:
main_2D_aniso_I: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue and with three occlusion sites as initial data, the results are shown in Figure 2.7.
main_2D_aniso_II: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue and with two occlusion sites as initial data, the results are shown in Figure 2.8.
main_2D_iso_I: Simulates the model (2.4.16 & 2.3.7) with isotropic tissue and with three occlusion sites as initial data, the results are shown in Figure 2.5.
main_2D_iso_II: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue and with two occlusion sites as initial data, the results are shown in Figure 2.6.
main_2D_lower_grades: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue, with two occlusion sites as initial data, and with stronger proton buffering the results are shown in Figure 2.9.
main_2D_modified_source: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue, with two occlusion sites as initial data, and with a modified source term for glioma cell equation the results are shown in Figure 2.10.



Method: for time discretisation IMEX method is used. Diffusion is solved implicitly(implicit Euler) while source and taxis/advection terms explicitly(Euler).
 
For space discretisation, Weikert's 3X3 (central diff, 9 points as mentioned in the Table 2.3 of the thesis) stencil discretisation(for anisotropic diffusion) has been used for tumor diffusion and standard 5 point stencil(again central diff, as shown in Table 2.2)  for acid diffusion. Upwind scheme(first order) is used in both x and y direction for the taxis term.

Time: 4-5 minutes(on 8 core, 16gb RAM, i7-processor @ 3.60GHz speed, Linux) to 5-6 minutes(dual core, 8GB RAM, i5-processor, @2,3GHz speed, Mac)

Requirement: Matlab(2018a or later version, it should work for older version also)
