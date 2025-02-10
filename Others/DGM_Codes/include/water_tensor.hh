#ifndef WATER_TENSOR_HH
#define WATER_TENSOR_HH

#include <deal.II/base/tensor_function.h>
#include <string>
#include <netcdf>
#include <array>
#include <vector>

using namespace dealii;
using namespace netCDF;
using namespace netCDF::exceptions;
template <int dim>
class WaterTensor : public dealii::TensorFunction<2,dim>
{
public:
  typedef dealii::TensorFunction<2,dim> BaseT;
  typedef typename BaseT::value_type value_type;
  WaterTensor(const std::string & filename);
  virtual value_type value(const Point<dim> & p)const;

private:
  std::array<double,3> m_spacing;
  std::array<unsigned int,5> m_sizes;
  std::vector<value_type> m_waterTensors;
};

template <int dim>
WaterTensor<dim>::WaterTensor(const std::string & filename):
  BaseT()
{
  std::cout << "Starting read DW"<<std::endl;
  NcFile file(filename,NcFile::read);
  NcVar spacing_var = file.getVar("spacing");
  spacing_var.getVar(&(m_spacing.front()));

  NcVar diff_var = file.getVar("DW");
  std::vector<NcDim> dims = diff_var.getDims();
  if(dims.size() != 5)
    throw std::runtime_error("Problem with dimensions");
  for(unsigned int i=0;i<5 ;++i)
    m_sizes[i] = dims.at(i).getSize();

  double * dt = new double[std::accumulate(m_sizes.begin(),
					   m_sizes.end(),
					   1,
					   std::multiplies<unsigned int>())];
  diff_var.getVar(dt);
  m_waterTensors.resize(m_sizes[0] * m_sizes[1] * m_sizes[2]);
  
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
	    (j + m_sizes[1] * i);
	  m_waterTensors.at(mapped_index) = tensor;
	}
  delete[] dt;
  std::cout << "Done reading DW" << std::endl;
}

template <int dim>
typename WaterTensor<dim>::value_type WaterTensor<dim>::
value(const Point<dim> & p)const
{
  const unsigned int slice = 100;
  std::array<unsigned int,3> index;
  index[2] = slice;
  for(unsigned int i=0;i<dim;++i)
    index[i] = static_cast<unsigned int>(p[i] / m_spacing[i] -
					 std::numeric_limits<float>::epsilon());

  return m_waterTensors.at(index[2] + m_sizes[2] * (index[1] + m_sizes[1] *index[0]));
}


#endif // WATER_TENSOR_HH
