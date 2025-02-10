//SGN
#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/lac/vector.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "biology_hypo_micro.hh"
#include "assembler_hypo_micro.hh"
using namespace dealii;
template <int dim>
class newS
{
public: 
	newS(const double dt, const unsigned int ref,
			const unsigned int fe_deg_);
	void run();
private:
	void setup();
	void solve_mass(const Vector<double> &rhs, Vector<double> &solution) const;
	void solve_diffusion_m(const Vector<double> &rhs, Vector<double> &solution) const;
	void solve_diffusion_s(const Vector<double> &rhs, Vector<double> &solution) const;
	void output_function(const Vector<double> &function,
			const std::string filename,
			const std::string varname) const;
	void advance_step(const std::vector<Vector<double> > &old_solution,
			std::vector<Vector<double> > &solution);
	Triangulation<dim> triangulation;
	FE_DGQ<dim> fe;
	DoFHandler<dim> dof_handler;

	Biology_s biology_s;
	Assembler_s<dim> assembler_s;
	const double time_step;
	const unsigned int refine;
	const unsigned int fe_deg;
	unsigned int dofs_per_cell;

	SparseMatrix<double> advance_diffusion_m;
	SparseMatrix<double> advance_diffusion_s;
	SparseMatrix<double> advance_advection;
	PreconditionJacobi<SparseMatrix<double> > diffusion_precond;

	Vector<double> rhs_m,rhs_s;
	std::vector<Vector<double> > source;
};
template <int dim>
newS<dim>::newS (const double dt,
		const unsigned int ref,
		const unsigned int fe_deg_):
		fe(fe_deg_),
		dof_handler(triangulation),
		time_step(dt),
		refine(ref),
		fe_deg(fe_deg_)
		{
		}


template <int dim>
void newS<dim>::setup()
{

	// reading the grid
	//	GridIn<dim> grid_in;
	//	grid_in.attach_triangulation (triangulation);
	//	std::ifstream input_file("grid.msh");
	//	grid_in.read_msh (input_file);

	Point<dim> p1;
	Point<dim> p2;
	p1(0) = 0;
	p1(1) = 0;
	p2(0) = 1000;
	p2(1) = 1000;

	GridGenerator::hyper_rectangle(triangulation, p1,p2);

	triangulation.refine_global(6);
	dof_handler.distribute_dofs(fe);
	dofs_per_cell = fe.dofs_per_cell;



	assembler_s.setup(dof_handler, biology_s);
	assembler_s.assemble_diffusion_matrix_m();
	assembler_s.assemble_diffusion_matrix_s();///

	//std::cout<<assembler_s.assemble_diffusion_matrix_m[0]<<std::endl;

	assembler_s.diffusion_matrix_m *= time_step;
	advance_diffusion_m.reinit(assembler_s.flux_sparsity);
	advance_diffusion_m.copy_from(assembler_s.mass_matrix);
	advance_diffusion_m.add(biology_s.D_m, assembler_s.diffusion_matrix_m);

	assembler_s.diffusion_matrix_s *= time_step;///
	advance_diffusion_s.reinit(assembler_s.flux_sparsity);
	advance_diffusion_s.copy_from(assembler_s.mass_matrix);
	advance_diffusion_s.add(biology_s.D_s, assembler_s.diffusion_matrix_s);


	advance_advection.reinit(assembler_s.flux_sparsity);
	rhs_m.reinit(dof_handler.n_dofs());
	rhs_s.reinit(rhs_m.size());
	source.resize(2);
	for (unsigned int i=0; i<2; ++i)
		source[i].reinit(rhs_m.size());
}//
template <int dim>
void newS<dim>::solve_mass(const Vector<double> &rhs,
		Vector<double> &solution) const
		{
	SolverControl solver_control(10000,1e-12);
	SolverCG<> solver(solver_control);
	solver.solve(assembler_s.mass_mass_matrix, solution, rhs, assembler_s.mass_precond);
		}

template <int dim>
void newS<dim>::solve_diffusion_m(const Vector<double> &rhs,
		Vector<double> &solution) const
		{
	SolverControl solver_control(10000,1e-12);
	SolverBicgstab<> solver(solver_control);
	solver.solve(advance_diffusion_m, solution, rhs, PreconditionIdentity());
		}
template <int dim>
void newS<dim>::solve_diffusion_s(const Vector<double> &rhs,
		Vector<double> &solution) const
		{
	SolverControl solver_control(10000,1e-12);
	SolverBicgstab<> solver(solver_control);
	solver.solve(advance_diffusion_s, solution, rhs, PreconditionIdentity());
		}

template <int dim>
void newS<dim>::output_function(const Vector<double> &function,
		const std::string filename,
		const std::string varname) const
		{
	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(function, varname);


	data_out.build_patches(assembler_s.fe_degree+2);// fe_degree+2 ??
	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
		}
template <int dim>
void newS<dim>::advance_step(const std::vector<Vector<double> > &old_solution,
		std::vector<Vector<double> > &solution)
		{
	advance_advection.copy_from(assembler_s.mass_matrix);
	//assembler_s.assemble_advection_matrix(old_solution[1],old_solution[1]);
	assembler_s.assemble_advection_matrix(old_solution[1]);

	assembler_s.advection_matrix *= -time_step;
	advance_advection.add(1.0, assembler_s.advection_matrix);

	advance_advection.vmult(rhs_m,old_solution[0]);
	assembler_s.mass_mass_matrix.vmult(rhs_s,old_solution[1]);


	assembler_s.assemble_source(old_solution,source);
	rhs_m.add(time_step,source[0]);
	rhs_s.add(time_step,source[1]);

	solve_diffusion_m(rhs_m, solution[0]);
	solve_diffusion_s(rhs_s, solution[1]);

		}
template <int dim>
void newS<dim>::run()
{
	Timer timer;
	std::cout<<"Setting up and assembling..."<<std::endl;
	timer.start();
	setup();

	std::cout<<"Assembly done."<<std::endl;
	timer.stop();
	std::cout << "Elapsed CPU time: " << timer() << " seconds."<<std::endl;
	std::cout << "Elapsed wall time: " << timer.wall_time() << " seconds."<<std::endl;
	std::cout << "Total computational cells: "
			<<triangulation.n_active_cells()
			<<std::endl;
	double Tend =25;
	std::cout<< "Advancing the solution..."<<std::endl;
	timer.reset();
	timer.start();

	std::vector<Vector<double> > solution(2);

	for (unsigned int i=0; i<2; ++i)
		solution[i].reinit(assembler_s.n_dofs);
	std::vector<Vector<double> > old_solution(solution);
	assembler_s.initial_values(old_solution);


	output_function(old_solution[0], "initial_M.vtk","M");
	output_function(old_solution[1], "initial_S.vtk","S");
	double time = time_step;
	for (unsigned int timestep_number=0; time<=Tend;
			time+=time_step, ++timestep_number)
	{
		advance_step(old_solution, solution);
		old_solution = solution;

		if (timestep_number%100==0)
		{
			output_function(solution[0], "out//m_solution"
					+Utilities::int_to_string(timestep_number,3)+".vtk",
					"M");
			output_function(solution[1], "out//s_solution"
					+Utilities::int_to_string(timestep_number,3)+".vtk",
					"S");
			//std::cout<< ".... Just reached  " << timestep_number << "  steps....."<<std::endl;


		}
		std::cout<<'\r'
						<<std::setw(5)<<time<<" /"
						<<std::setw(5)<<Tend<<std::flush;

	}

	timer.stop();
	std::cout<< "Simulation done"<<std::endl;
	std::cout << "Elapsed CPU time: " << timer() << " seconds."<<std::endl;
	std::cout << "Elapsed wall time: " << timer.wall_time() << " seconds."<<std::endl;
}

int main ()
{
	try
	{
		using namespace dealii;

		deallog.depth_console (0);

		const double dt = 1e-2;
		newS<2> main(dt,0,1);
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
