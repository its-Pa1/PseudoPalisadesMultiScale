# Difference Directory

## Main Files/Scripts

- **main_diff.m**: The main MATLAB script which calls `main_2D.m` and `main_2D_without_g.m`. These two scripts simulate the model (2.4.16 & 2.3.7) with and without g(S) respectively. The produced results are shown in Figure 2.11 of the thesis.

## Functions

- **tissue_Q_macro.m**: Computes the macroscopic tissue Q described in subsection 2.5.2 of the thesis.
- **function_g.m**: Computes the g(S(x)) involved in equation 2.4.16.
- **tissue_q_un.m**: Computes the undirected tissue small q given by equation 2.5.4.
- **set_tumor_diff.m**: Assembles the tumor diffusion matrix and stores other related variables.
- **set_acidity_diff.m**: Assembles the tumor diffusion matrix.

## Method

- **Time Discretization**: IMEX method is used. Diffusion is solved implicitly (implicit Euler) while source and taxis/advection terms are solved explicitly (Euler).
- **Space Discretization**: 
  - Weikert's 3x3 stencil discretization (central difference, 9 points as mentioned in Table 2.3 of the thesis) is used for tumor diffusion.
  - Standard 5-point stencil (central difference, as shown in Table 2.2) is used for acid diffusion.
  - Upwind scheme (first order) is used in both x and y directions for the taxis term (see subsection 2.5.3 of the thesis).

## Performance

- **Time**: 
  - 8-10 minutes on an 8-core, 16GB RAM, i7-processor @ 3.60GHz speed, Linux.
  - 10-12 minutes on a dual-core, 8GB RAM, i5-processor @ 2.3GHz speed, Mac.

## Requirements

- **MATLAB**: 2018a or later version (it should work for older versions as well).