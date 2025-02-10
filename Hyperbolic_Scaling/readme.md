# Hyperbolic_Scaling Directory

## Main Files/Scripts

- **main_hyp_all.m**: The main MATLAB script which calls all the functions (scripts) and simulates the 2-D model on a square grid. It produces most of the results shown in subsection 2.5.5 of "Diss_Kumar_Pawan.pdf" present in the parent directory.

## Functions

- **tissue_Q_macro.m**: Computes the macroscopic tissue Q described in subsection 2.5.2 of the thesis.
- **function_g.m**: Computes the g(S(x)) involved in equation 2.4.46.
- **tissue_unimodal.m**: Computes the directed tissue (q) given by equation 2.5.5.
- **set_tumor_diff.m**: Assembles the tumor diffusion matrix and stores other related variables.
- **set_acidity_diff.m**: Assembles the tumor diffusion matrix.

## Scripts

- **main_SU_delta1.m**: Simulates the model (2.4.46 & 2.3.7) with anisotropic tissue and with two occlusion sites as initial data and with delta = 1 in equation 2.5.5. The results are shown in Figures 2.12 and 2.13.
- **main_SU_exp2_delta0p2.m**: Simulates the model (2.4.46 & 2.3.7) with anisotropic tissue and with two occlusion sites as initial data and with delta = 0.2 in equation 2.5.5. The results are shown in Figures 2.14 and 2.15.

## Method

- **Time Discretization**: IMEX method is used. Diffusion is solved implicitly (implicit Euler) while source and taxis/advection terms are solved explicitly (Euler).
- **Space Discretization**: 
  - Weikert's 3x3 stencil discretization (central difference, 9 points as mentioned in Table 2.3 of the thesis) is used for tumor diffusion.
  - Standard 5-point stencil (central difference, as shown in Table 2.2) is used for acid diffusion.
  - Second order upwind scheme with Van Leer flux limiter is used in both x and y directions for the taxis term (see subsection 2.5.3 of the thesis).

## Performance

- **Time**: 
  - 30-35 minutes on an 8-core, 16GB RAM, i7-processor @ 3.60GHz speed, Linux.
  - 35-40 minutes on a dual-core, 8GB RAM, i5-processor @ 2.3GHz speed, Mac.

## Requirements

- **MATLAB**: 2018a or later version (it should work for older versions as well).