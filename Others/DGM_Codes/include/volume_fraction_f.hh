#ifndef VOLUME_FRACTION_HH
#define VOLUME_FRACTION_HH


#include <deal.II/base/function.h>
#include <string>
#include <netcdf>
#include <array>
#include <vector>
#include <numeric>
using namespace dealii;
using namespace netCDF;
using namespace netCDF::exceptions;
template <int dim>
class VolumeFraction :public dealii::Function<dim>
{
public:
  typedef dealii::Function<dim> BaseT;
  //typedef typename BaseT;
  VolumeFraction(const std::string & filename);
  virtual double value(const Point<dim> & p,
		       const unsigned int  component = 0)const;

private:
  std::array<double,3> m_spacing;
  std::array<unsigned int,3> m_sizes;
  std::vector<double> m_volume_fraction;
};

template <int dim>
VolumeFraction<dim>::VolumeFraction(const std::string & filename):BaseT()
{
  std::cout << "Starting read Q"<<std::endl;
  NcFile file(filename,NcFile::read);
  NcVar spacing_var = file.getVar("spacing");
  spacing_var.getVar(&(m_spacing.front()));

  NcVar Q_var = file.getVar("Q");
  std::vector<NcDim> dims = Q_var.getDims();
  if(dims.size() != 3)
    throw std::runtime_error("Problem with sizes of Q");
  for(unsigned int i=0;i<3 ;++i)
    m_sizes[i] = dims.at(i).getSize();

  double * q = new double[std::accumulate(m_sizes.begin(),
					   m_sizes.end(),
					   1,
					   std::multiplies<unsigned int>())];
  Q_var.getVar(q);
  m_volume_fraction.resize(m_sizes[0] * m_sizes[1] * m_sizes[2]);
  
  for(unsigned int i=0;i<m_sizes[0];++i)
    for(unsigned int j=0;j<m_sizes[1];++j)
      for(unsigned int k=0;k<m_sizes[2];++k)
	{
	  
	  const unsigned int mapped_index =
	    k + m_sizes[2] *
	    (j + m_sizes[1] * i);
	  m_volume_fraction.at(mapped_index) = q[mapped_index];
	}
  delete[] q;
  std::cout << "Done reading Q" << std::endl;
}

template <int dim>
double VolumeFraction<dim>::value(const Point<dim> & p,
				   const unsigned int /* component*/)const
{
  const unsigned int slice = 100;
  std::array<unsigned int,3> index;
  index[2] = slice;
  for(unsigned int i=0;i<dim;++i)
    index[i] = static_cast<unsigned int>(p[i] / m_spacing[i] -
					 std::numeric_limits<float>::epsilon());

  return m_volume_fraction.at(index[2] + m_sizes[2] * (index[1] + m_sizes[1] *index[0]));
}


#endif // VOLUME_FRACTION_HH
