#ifndef WATER_TENSOR_HH
#define WATER_TENSOR_HH
// This header file deals with the water diffusion tensor obtained from DTI data.
// all indeing are row-major

//To use the tensor library from deal.II
#include <deal.II/base/tensor_function.h>

// other libraries from std and netcdf
#include <string>
#include <netcdf>
#include <array>
#include <vector>

// declaration of namespace to avoid calling them like deal.II::
using namespace dealii;
using namespace netCDF;
using namespace netCDF::exceptions;

// water diffusion tensor class
template <int dim>
class WaterTensor : public dealii::TensorFunction<2,dim> // dim dimensional tensor
{
public:
	typedef dealii::TensorFunction<2,dim> BaseT; //
	typedef typename BaseT::value_type value_type;
	WaterTensor(const std::string & filename); // file destination
	virtual value_type value(const Point<dim> & p)const; // virtual function to calculate DW at called point

private:
	std::array<double,3> m_spacing; // size of each voxels here its 2 in each direction
	std::array<unsigned int,5> m_sizes; // XxYxZx3x3 or XxYxZxdimxdim
	std::vector<value_type> m_waterTensors; // array of water tensor
};

template <int dim>
WaterTensor<dim>::WaterTensor(const std::string & filename):
BaseT()
{

	// netcdf stuffs: reading the data
	std::cout << "Reading DW"<<std::endl;
	NcFile file(filename,NcFile::read);
	NcVar spacing_var = file.getVar("spacing");
	spacing_var.getVar(&(m_spacing.front()));
	NcVar diff_var = file.getVar("DW");

	// check whether the dimension is correct
	std::vector<NcDim> dims = diff_var.getDims();
	if(dims.size() != 5)
		throw std::runtime_error("Problem with dimensions");

	//obtain the size of data
	for(unsigned int i=0;i<5 ;++i)
		m_sizes[i] = dims.at(i).getSize();

	// allocate the memory XxYxZxdimxdim and arrange the data as an 1D array
	double * dt = new double[std::accumulate(m_sizes.begin(),
			m_sizes.end(),
			1,
			std::multiplies<unsigned int>())];
	diff_var.getVar(dt);

	m_waterTensors.resize(m_sizes[0] * m_sizes[1] * m_sizes[2]);// memory for one dimensional array of Dw at each point in XYZ space

	//linear indexing in 5-D
	for(unsigned int i=0;i<m_sizes[0];++i)
		for(unsigned int j=0;j<m_sizes[1];++j)
			for(unsigned int k=0;k<m_sizes[2];++k)
			{
				value_type tensor;
				for(unsigned int ii=0;ii<dim;++ii)
					for(unsigned int jj=0;jj<dim;++jj)
					{
						const unsigned int index =
								jj + m_sizes[4] *
								( ii + m_sizes[3] *
										( k + m_sizes[2] *
												(j + m_sizes[1] * i)));
						tensor[ii][jj] = dt[index];
					}
				const unsigned int mapped_index =
						k + m_sizes[2] *
						(j + m_sizes[1] * i);// again linear index in 3-d

				m_waterTensors.at(mapped_index) = tensor; // store the tensor at each index
			}
	delete[] dt; // free the allocated memory
	std::cout << "Done reading DW" << std::endl;
}
// this function calculate the water tensor at a called point
// it calls the stored tensor at the indexed point considering the spacing
template <int dim>
typename WaterTensor<dim>::value_type WaterTensor<dim>::
value(const Point<dim> & p)const
{
	const unsigned int slice = 100; // for 2-D tensors, it fix the Z cordinate
	std::array<unsigned int,3> index;
	index[2] = slice;
	// in case of 3-D, the loop overwrites the assigned 3rd enrty
	for(unsigned int i=0;i<dim;++i){
		index[i] = static_cast<unsigned int>(p[i] / m_spacing[i] +
				std::numeric_limits<float>::epsilon());

	}
	// finally it desired water tensor at point p
	return m_waterTensors.at(index[2] + m_sizes[2] * (index[1] + m_sizes[1] *index[0]));
}


#endif // WATER_TENSOR_HH
