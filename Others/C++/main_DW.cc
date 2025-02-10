// This is the main file, which calls the header file and computes
// Dw at called points


// some useful libraries and the DW header file
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include "water_tensor.hh"


#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace dealii;

template<int dim>
class DTI_data
{
public:

	DTI_data(Point<dim> p); //class constructor
	void run();

private:
	Point<dim> pt ; // point at which Dw is required
	typedef WaterTensor<dim> Water_Tensor; // water tensor object
	Water_Tensor water_tensor;
};
// reding the data and space-point pt
template<int dim>
DTI_data<dim>::DTI_data(Point<dim> p)
:pt(p),
 water_tensor("data/data.ncdf")
 {}
// this function call the value function and output the Dw value at the point
template<int dim>
void DTI_data<dim>::run(){
	Tensor<2,dim> DW;
	DW = water_tensor.value(pt);
	std::cout<<"The W value at [" << pt << "] = "<<DW<< std::endl;
}

int main ()
{
	try
	{
		using namespace dealii;

		deallog.depth_console (0); // deal.II stuff
		// 3D point
		Point<3> p ;
		p(0) = 120.0;
		p(1) = 180.0;
		p(2) = 200.0;
		DTI_data<3> main(p);


		//2D point
		//Point<2> p ;
		//p(0) = 100.0;
		//p(1) = 150.0;
		//DTI_data<2> main(p);

		main.run();
	}
	catch (std::exception &exc)
	{
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
				<< exc.what() << std::endl
				<< "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;

		return 1;
	}
	catch (...)
	{
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
				<< "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		return 1;
	}

	return 0;
}
