# DGM_Codes Directory

## Overview

This directory contains the C++ codes for the Discontinuous Galerkin Method (DGM) used in the multiscale modeling of glioma pseudopalisades. The codes are designed to solve partial differential equations (PDEs) related to tumor growth and diffusion processes.

## Main Files/Scripts

- **main.cc**: The main C++ file which calls all the functions and simulates the model. It includes the necessary headers and sets up the simulation environment.

- **assembler_hypo_micro.hh**: This header file defines the `Assembler_s` class, which is responsible for assembling various matrices (diffusion, advection, and source) and computing initial values for the simulation. It includes methods for setting up the problem, assembling matrices, and solving the mass matrix.

- **diffusion_tensor.hh**: This header file defines the `DiffusionTensor` class, which reads the diffusion tensor data from a NetCDF file and provides the tensor values at specified points.

- **volume_fraction_f.hh**: This header file defines the `VolumeFraction` class, which reads the volume fraction data from a NetCDF file and provides the volume fraction values at specified points.

- **water_tensor.hh**: This header file defines the `WaterTensor` class, which reads the water diffusion tensor data from a NetCDF file and provides the tensor values at specified points.

- **biology_hypo_micro.hh**: This header file defines the `Biology_s` class, which includes functions for computing source terms and initial values for the biological model.

## Functions

- **setup**: Initializes the degrees of freedom (DoF) handler, biology parameters, and sparsity patterns for the matrices.
- **assemble_diffusion_matrix_m**: Assembles the diffusion matrix for the model.
- **assemble_diffusion_matrix_s**: Assembles the diffusion matrix for the source term.
- **assemble_source**: Assembles the source term matrix.
- **assemble_advection_matrix**: Assembles the advection matrix.
- **initial_values**: Computes the initial values for the simulation.
- **solve_mass**: Solves the mass matrix equation.

## Method

- **Space Discretization**: Symmetric Interior Penalty Galerkin (SIPG) method is used.
- **Time Discretization**: IMEX (Implicit-Explicit) scheme is used.

## Important Note

The user has to keep the `data.ncdf` file in the `data` directory. The author of this repository doesn't have the copyright of the data to put here.

## Required Libraries

- deal.II (version 9.1.0 or later)
- netCDF4
- HDF5

## Installation

1. **Install deal.II**: Follow the instructions on the [deal.II website](https://www.dealii.org/) to install deal.II.
2. **Link deal.II using CMake**: Ensure that deal.II is properly linked using CMake. You can use the following CMake configuration:

    ```bash
    cmake -DDEAL_II_DIR=/path/to/deal.II .
    ```

## Usage

1. Ensure that the required libraries are installed.
2. Place the `data.ncdf` file in the `data` directory.
3. Compile and run the main C++ program `hypo_macro.cc`.

## Concept

The idea behind this code is to use the Discontinuous Galerkin Method (DGM) to solve PDEs related to tumor growth and diffusion processes. The code reads the diffusion tensor data from a NetCDF file and uses it to assemble the necessary matrices for the simulation.

## Important Considerations

One should be careful about the alignment/synchronization of the brain slice (grid) and the data. This data fits well with the grid (file: `grid.msh`) displayed in the file `brain_slice.jpg`.

One of the benefits of using deal.II/DUNE libraries is that they automatically detect the required computational domain (also degree of freedom, etc.) from the grid. Of course, this can be done without those libraries as well.
