# C++ Directory

## Data

The DTI data (`C++/data/data.ncdf`) is stored in a NetCDF file.

### Viewing/Opening the Data

To view/open/use this file in a program (C++), the following libraries are required:
- netCDF4 (or later version)
- HDF5
- HDFView

This file can also be viewed/opened using MATLAB:
- To know about the data: `ncdisp('data.ncdf');`
- To load the data: `ncread('data.ncdf', 'DW');`

The `ncdisp` function shows everything about the data in the command window.

### Data Structure

The DW data is stored in the following way: `XxYxZx3x3` (but MATLAB reads it the other way around).

### Important Note

The user has to keep the `data.ncdf` file in the `data` directory. The author of this repository doesn't have the copyright of the data to copy here.

## C++ Program

### Required Libraries

- CMake (version 2.8.8 or later)
- deal.II (version 9.1.0 or later)
- netCDF4
- HDF5

### Main File

The main file `C++/main_DW.cc` calls the header file `C++/water_tensor.hh` and calculates the water diffusion tensor at a specified point.

### Concept

The idea behind this code is to store the entire DW data as an array and then evaluate it in the desired dimension and at the desired point using linear indexing.

### Important Considerations

One should be careful about the alignment/synchronization of the brain slice (grid) and this data. This data fits well with the grid (file: `grid.msh`) displayed in the file `brain_slice.jpg`.

One of the benefits of using deal.II/DUNE libraries is that they automatically detect the required computational domain (also degree of freedom, etc.) from the grid. Of course, this can be done without those libraries as well.