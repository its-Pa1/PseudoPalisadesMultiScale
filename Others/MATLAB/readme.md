# MATLAB Directory

## Main Files/Scripts

- **WaterTensor.m**: Reads and stores the DW data from the `data.ncdf` file and plots the fractional anisotropy of DW for a given domain.

First, the DW data is loaded and then converted from 3x3xZxYxX to XxYxZx3x3. After that, the entire 5D matrix is stored in a one-dimensional array using row-major linear indexing, as MATLAB stores it as column-major indexing. Then DW at each point in the domain is stored (depending upon the chosen dimension 2 or 3), followed by the FA calculation and its plotting.

- **WaterTensor2.m**: Performs the same task as `WaterTensor.m`. In this script, the data has not been transformed and is used the way it is stored by MATLAB's `ncread` function. The MATLAB-defined function `sub2ind` is used to work with column-major indexing.

- **WaterTensor3.m**: Also performs the same task. It is the simplest and fastest script. It doesn't contain any indexing. It just accesses the DW for a given domain and then trims the 3D data to 2D.

## Requirements

- **MATLAB**: 2020a/2018a (should work on older versions as well).

### Important Note

The user has to keep the `data.ncdf` file in this` directory. The author of this repository doesn't have the copyright of the data to put here.