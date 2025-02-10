# Parabolic Directory

## Main Files/Scripts

- **main_all.m**: The main MATLAB script which calls all the functions and simulates the 2-D model on a square grid. It produces most of the results shown in subsection 2.5.4 of "Diss_Kumar_Pawan.pdf" present in the parent directory.

## Functions

- **tissue_Q_macro.m**: Computes the macroscopic tissue Q described in subsection 2.5.2 of the thesis.
- **function_g.m**: Computes the g(S(x)) involved in equation 2.4.16.
- **tissue_q_un.m**: Computes the undirected tissue small q given by equation 2.5.4.
- **set_tumor_diff.m**: Assembles the tumor diffusion matrix and stores other related variables.
- **set_acidity_diff.m**: Assembles the tumor diffusion matrix.

## Scripts

- **main_2D_aniso_I.m**: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue and with three occlusion sites as initial data. The results are shown in Figure 2.7.
- **main_2D_aniso_II.m**: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue and with two occlusion sites as initial data. The results are shown in Figure 2.8.
- **main_2D_iso_I.m**: Simulates the model (2.4.16 & 2.3.7) with isotropic tissue and with three occlusion sites as initial data. The results are shown in Figure 2.5.
- **main_2D_iso_II.m**: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue and with two occlusion sites as initial data. The results are shown in Figure 2.6.
- **main_2D_lower_grades.m**: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue, with two occlusion sites as initial data, and with stronger proton buffering. The results are shown in Figure 2.9.
- **main_2D_modified_source.m**: Simulates the model (2.4.16 & 2.3.7) with anisotropic tissue, with two occlusion sites as initial data, and with a modified source term for the glioma cell equation. The results are shown in Figure 2.10.

## Method

- **Time Discretization**: IMEX method is used. Diffusion is solved implicitly (implicit Euler) while source and taxis/advection terms are solved explicitly (Euler).
- **Space Discretization**: 
  - Weikert's 3x3 stencil discretization (central difference, 9 points as mentioned in Table 2.3 of the thesis) is used for tumor diffusion.
  - Standard 5-point stencil (central difference, as shown in Table 2.2) is used for acid diffusion.
  - Upwind scheme (first order) is used in both x and y directions for the taxis term.

## Performance

- **Time**: 
  - 4-5 minutes on an 8-core, 16GB RAM, i7-processor @ 3.60GHz speed, Linux.
  - 5-6 minutes on a dual-core, 8GB RAM, i5-processor @ 2.3GHz speed, Mac.

## Requirements

- **MATLAB**: 2018a or later version (it should work for older versions as well).