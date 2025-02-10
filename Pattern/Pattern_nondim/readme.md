# Pattern_nondim Directory

## Main Files/Scripts

- **Pattern_main_const.m**: The main MATLAB script which calls all the functions and simulates the non-dimensional 1-D model with constant diffusion coefficient D. It plots tumor and acidity with respect to time and space.
- **Pattern_main_D1.m**: The main MATLAB script which calls all the functions and simulates the non-dimensional 1-D model with diffusion coefficient D1 as mentioned in subsection 2.5.6 of the thesis. It plots tumor and acidity with respect to time and space.
- **Pattern_main_D2.m**: The main MATLAB script which calls all the functions and simulates the non-dimensional 1-D model with diffusion coefficient D2 as mentioned in subsection 2.5.6 of the thesis. It plots tumor and acidity with respect to time and space.

## Functions

- **set_diff_mat_const.m**: Assembles the tumor diffusion matrix and stores constant D.
- **set_diff_mat_D1.m**: Assembles the tumor diffusion matrix and stores D1.
- **set_diff_mat_D2.m**: Assembles the tumor diffusion matrix and stores D2.

## Method

- **Time Discretization**: IMEX method is used. Diffusion is solved implicitly (implicit Euler) while source and taxis/advection terms are solved explicitly (Euler).
- **Space Discretization**: 
  - Standard 5-point stencil (central difference) is used for diffusion.
  - Upwind scheme (first order) is used for the taxis term.

## Requirements

- **MATLAB**: 2018a or later version (it should work for older versions as well).