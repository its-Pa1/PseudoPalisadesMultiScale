//SGN
#ifndef ASSEMBLER_HYPO_MICRO_HH
#define ASSEMBLER_HYPO_MICRO_HH

#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "biology_hypo_micro.hh"
#include "include/diffusion_tensor.hh"
#include "include/volume_fraction_f.hh"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace dealii;
template <int dim>
class Assembler_s
{
public:
	Assembler_s();
	void setup(const DoFHandler<dim> &dof_handler_,
			const Biology_s &biology_s_);
	void assemble_diffusion_matrix_m();
	void assemble_diffusion_matrix_s();
	void assemble_source(const std::vector<Vector<double> > &functions,
			std::vector<Vector<double> > &source);
	void assemble_advection_matrix(const Vector<double> &h);
	void initial_values(std::vector<Vector<double> > &initial) const;

	DoFHandler<dim> const *dof_handler;
	FiniteElement<dim> const *fe;
	Biology_s const *biology_s;
	unsigned int n_dofs;
	unsigned int dofs_per_cell;
	unsigned int fe_degree;
	SparsityPattern flux_sparsity;
	SparsityPattern mass_sparsity;
	SparseMatrix<double> diffusion_matrix_m;
	SparseMatrix<double> diffusion_matrix_s;
	SparseMatrix<double> mass_matrix;
	SparseMatrix<double> mass_mass_matrix;
	SparseMatrix<double> advection_matrix;
	PreconditionJacobi<SparseMatrix<double> > mass_precond;
private:
	////////////////////////////////////////////////////////////////////////////////////////////////////
	void compute_p_d(const std::vector<double> &points,
			const std::vector<double> &func_values,
			const std::vector<bool> &sides,
			std::vector<double>  &derivatives) const;
	void compute_drift_field(const std::vector<Point<dim> > &points,
			const std::vector<Tensor<2,dim> > &DT_values,
			const std::vector<bool> &sides,
			std::vector<Tensor<1,dim> > &drifts) const;
	void compute_sub_func(const double &Qx,
			double &g);
	typedef DiffusionTensor<2> DT;
	DT diffusionTensor;

	//	typedef VolumeFraction<2> Q;
	//	Q volumeFraction;

	const double scale_hour = 60*60;
	const double scale_day = scale_hour*24;
	const double scale_week = scale_day*7;
	const double scale_month = scale_day*30;
	const double speed = 50.0/3600;

	const double scale = scale_month; // 1 = per sec
	const double s = speed*scale;//50mu/hr
	const double lambda_0 = 0.8*scale;
	const double lambda_1 = 0.1*scale;
	const double k_new = 0.1;
	const double T_scale = pow(s,2)/lambda_0;

//	const double R_0 = 1e+5;
//	const double k_p = 0.1*scale;
//	const double k_m = 0.1*scale;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	struct Scratch_diff
	{
		Scratch_diff(const FiniteElement<dim> &fe,
				const Quadrature<dim> &quadrature,
				const Quadrature<dim-1> &face_quadrature);
		Scratch_diff(const Scratch_diff &scratch);

		FEValues<dim> fe_values;
		FEFaceValues<dim> fe_face_values;
		FEFaceValues<dim> fe_face_values_neighbor;

		FullMatrix<double> u1_v1;
		FullMatrix<double> u2_v2;
		FullMatrix<double> u1_v2;
		FullMatrix<double> u2_v1;
	};
	struct Copy_diff
	{
		FullMatrix<double> local_matrix;
		std::vector<std::vector<FullMatrix<double> > > neighbor_u_v;
		std::vector<std::vector<types::global_dof_index> > neighbor_dof_indices;
		std::vector<types::global_dof_index> local_dof_indices;
		std::vector<bool> sides;
	};
	void assemble_cell_diff(const typename DoFHandler<dim>::active_cell_iterator &cell,
			Scratch_diff &scratch, Copy_diff &copy);
	void copy_cell_diff(const Copy_diff &copy);

	//////////////////////////////////////////////////////////////////////////
	struct Scratch_diff2
	{
		Scratch_diff2(const FiniteElement<dim> &fe,
				const Quadrature<dim> &quadrature,
				const Quadrature<dim-1> &face_quadrature);
		Scratch_diff2(const Scratch_diff2 &scratch);

		FEValues<dim> fe_values;
		FEFaceValues<dim> fe_face_values;
		FEFaceValues<dim> fe_face_values_neighbor;

		FullMatrix<double> u1_v1;
		FullMatrix<double> u2_v2;
		FullMatrix<double> u1_v2;
		FullMatrix<double> u2_v1;
	};
	struct Copy_diff2
	{
		FullMatrix<double> local_matrix;
		std::vector<std::vector<FullMatrix<double> > > neighbor_u_v;
		std::vector<std::vector<types::global_dof_index> > neighbor_dof_indices;
		std::vector<types::global_dof_index> local_dof_indices;
		std::vector<bool> sides;
	};
	void assemble_cell_diff2(const typename DoFHandler<dim>::active_cell_iterator &cell,
			Scratch_diff2 &scratch, Copy_diff2 &copy);
	void copy_cell_diff2(const Copy_diff2 &copy);
	//////////////////////////////////////////////////////////////////////////
	struct Scratch_source
	{
		Scratch_source(const FiniteElement<dim> &fe,
				const Quadrature<dim> &quadrature);
		Scratch_source(const Scratch_source &scratch);

		FEValues<dim> fe_values;
		std::vector<std::vector<double> > function_values;
	};
	struct Copy_source
	{
		std::vector<Vector<double> > source_values;
		std::vector<types::global_dof_index> local_dof_indices;
	};
	void assemble_cell_source(const typename DoFHandler<dim>::active_cell_iterator &cell,
			Scratch_source &scratch, Copy_source &copy,
			const std::vector<Vector<double> > &functions);
	void copy_cell_source(const Copy_source &copy,
			std::vector<Vector<double> > &source);
	struct Scratch_adv
	{
		Scratch_adv(const FiniteElement<dim> &fe,
				const Quadrature<dim> &quadrature,
				const Quadrature<dim-1> &face_quadrature);
		Scratch_adv(const Scratch_adv &scratch);

		FEValues<dim> fe_values;
		FEValues<dim> fe_values_neighbor;
		FEFaceValues<dim> fe_face_values;
		FEFaceValues<dim> fe_face_values_neighbor;

		FullMatrix<double> u1_v1;
		FullMatrix<double> u2_v2;
		FullMatrix<double> u1_v2;
		FullMatrix<double> u2_v1;
	};
	struct Copy_adv
	{
		FullMatrix<double> local_matrix;
		std::vector<std::vector<FullMatrix<double> > > neighbor_u_v;
		std::vector<std::vector<types::global_dof_index> > neighbor_dof_indices;
		std::vector<types::global_dof_index> local_dof_indices;
		std::vector<bool> sides;
	};
	void assemble_cell_adv(const typename DoFHandler<dim>::active_cell_iterator &cell,
			Scratch_adv &scratch, Copy_adv &copy,
			const Vector<double> &h);
	//,const Vector<double> &v);
	void copy_cell_adv(const Copy_adv &copy);


	void solve_mass(const Vector<double> &rhs, Vector<double> &solution) const;
};

template <int dim>
Assembler_s<dim>::Assembler_s():
diffusionTensor("data/data.ncdf","DT")
//,
//volumeFraction("data/data.ncdf")
{}
template <int dim>
void Assembler_s<dim>::setup(const DoFHandler<dim> &dof_handler_,
		const Biology_s &biology_s_)
		{
	dof_handler = &dof_handler_;
	biology_s = &biology_s_;
	fe = &dof_handler_.get_fe();
	n_dofs = dof_handler->n_dofs();
	dofs_per_cell = fe->dofs_per_cell;
	fe_degree = fe->degree;

	DynamicSparsityPattern dsp(n_dofs);
	DoFTools::make_flux_sparsity_pattern (*dof_handler, dsp);
	flux_sparsity.copy_from (dsp);

	DynamicSparsityPattern dsp_m(n_dofs);
	DoFTools::make_sparsity_pattern(*dof_handler,dsp_m);
	mass_sparsity.copy_from(dsp_m);

	diffusion_matrix_m.reinit(flux_sparsity);
	diffusion_matrix_s.reinit(flux_sparsity);

	mass_matrix.reinit(flux_sparsity);
	mass_mass_matrix.reinit(mass_sparsity);
	advection_matrix.reinit(flux_sparsity);

	MatrixCreator::create_mass_matrix(*dof_handler,
			QGauss<dim>(fe_degree+1),
			mass_mass_matrix);
	MatrixCreator::create_mass_matrix(*dof_handler,
			QGauss<dim>(fe_degree+1),
			mass_matrix);
	mass_precond.initialize(mass_mass_matrix,1.0);
		}

template <int dim>
Assembler_s<dim>::Scratch_diff::Scratch_diff(const FiniteElement<dim> &fe,
		const Quadrature<dim> &quadrature,
		const Quadrature<dim-1> &face_quadrature)
:
		fe_values(fe, quadrature,
				update_values | update_gradients|
				update_quadrature_points| update_JxW_values),
				fe_face_values(fe, face_quadrature,
						update_values| update_gradients| update_normal_vectors|
						update_quadrature_points | update_JxW_values),
						fe_face_values_neighbor(fe, face_quadrature,
								update_values| update_gradients| update_normal_vectors|
								update_quadrature_points | update_JxW_values),
								u1_v1(fe.dofs_per_cell, fe.dofs_per_cell),
								u2_v2(fe.dofs_per_cell, fe.dofs_per_cell),
								u1_v2(fe.dofs_per_cell, fe.dofs_per_cell),
								u2_v1(fe.dofs_per_cell, fe.dofs_per_cell)
{}

template <int dim>
Assembler_s<dim>::Scratch_diff::Scratch_diff(const Scratch_diff &scratch)
:
fe_values(scratch.fe_values.get_fe(),
		scratch.fe_values.get_quadrature(),
		scratch.fe_values.get_update_flags()),
		fe_face_values(scratch.fe_face_values.get_fe(),
				scratch.fe_face_values.get_quadrature(),
				scratch.fe_face_values.get_update_flags()),
				fe_face_values_neighbor(scratch.fe_face_values_neighbor.get_fe(),
						scratch.fe_face_values_neighbor.get_quadrature(),
						scratch.fe_face_values_neighbor.get_update_flags()),
						u1_v1(scratch.u1_v1),
						u2_v2(scratch.u2_v2),
						u1_v2(scratch.u1_v2),
						u2_v1(scratch.u2_v1)
{}

template <int dim>
void Assembler_s<dim>::assemble_cell_diff(const typename DoFHandler<dim>::active_cell_iterator &cell,
		Scratch_diff &scratch, Copy_diff &copy)
		{
	const unsigned int n_q_points = scratch.fe_values.get_quadrature().size();
	const unsigned int n_face_q_points = scratch.fe_face_values.get_quadrature().size();

	copy.local_matrix.reinit(dofs_per_cell, dofs_per_cell);
	copy.neighbor_u_v.resize(4);
	copy.neighbor_dof_indices.resize(4);
	copy.sides.resize(4);
	copy.local_dof_indices.resize(dofs_per_cell);
	cell->get_dof_indices(copy.local_dof_indices);
	//const double edge_length = sqrt(0.5)*cell->diameter();

	scratch.fe_values.reinit(cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
			for (unsigned int q=0; q<n_q_points; ++q)
				copy.local_matrix(i,j) += scratch.fe_values.shape_grad(i,q)*
				scratch.fe_values.shape_grad(j,q)*scratch.fe_values.JxW(q);

	for (unsigned int face_n=0; face_n<4; ++face_n)
	{
		scratch.u1_v1 = 0;
		scratch.u1_v2 = 0;
		scratch.u2_v1 = 0;
		scratch.u2_v2 = 0;

		copy.neighbor_dof_indices[face_n].resize(dofs_per_cell);
		copy.neighbor_u_v[face_n].resize(4);
		for (unsigned int i=0; i<4; ++i)
			copy.neighbor_u_v[face_n][i].reinit(dofs_per_cell, dofs_per_cell);
		copy.sides[face_n] = false;

		if (!cell->at_boundary(face_n))
		{
			scratch.fe_face_values.reinit(cell,face_n);
			typename DoFHandler<dim>::active_cell_iterator
			neighbor = cell->neighbor(face_n);
			unsigned int neighbor_face = cell->neighbor_of_neighbor(face_n);
			scratch.fe_face_values_neighbor.reinit(neighbor,neighbor_face);
			copy.sides[face_n] = true;
			const double penalty = 100;//PENALTY HERE
			for (unsigned int q=0; q<n_face_q_points; ++q)
			{
				const Tensor<1,dim> normal_v = scratch.fe_face_values.normal_vector(q);
				for (unsigned int i=0; i<dofs_per_cell; ++i)
				{
					const Tensor<1,dim> phi_grad_i = scratch.fe_face_values.shape_grad(i,q);
					const double phi_i = scratch.fe_face_values.shape_value(i,q);
					const Tensor<1,dim> phi_grad_i_n = scratch.fe_face_values_neighbor.shape_grad(i,q);
					const double phi_i_n = scratch.fe_face_values_neighbor.shape_value(i,q);
					for (unsigned int j=0; j<dofs_per_cell; ++j)
					{
						const Tensor<1,dim> phi_grad_j = scratch.fe_face_values.shape_grad(j,q);
						const double phi_j = scratch.fe_face_values.shape_value(j,q);
						const Tensor<1,dim> phi_grad_j_n = scratch.fe_face_values_neighbor.shape_grad(j,q);
						const double phi_j_n = scratch.fe_face_values_neighbor.shape_value(j,q);

						scratch.u1_v1(i,j) += -0.5*(phi_grad_i*phi_j+phi_grad_j*phi_i)*normal_v*
								scratch.fe_face_values.JxW(q)+
								penalty*(phi_i*phi_j)*scratch.fe_face_values.JxW(q);
						scratch.u1_v2(i,j) += -0.5*(-phi_grad_j*phi_i_n+phi_grad_i_n*phi_j)*normal_v*
								scratch.fe_face_values.JxW(q)-
								penalty*(phi_i_n*phi_j)*scratch.fe_face_values.JxW(q);
						scratch.u2_v1(i,j) += -0.5*(phi_grad_j_n*phi_i-phi_grad_i*phi_j_n)*normal_v*
								scratch.fe_face_values.JxW(q)-
								penalty*(phi_i*phi_j_n)*scratch.fe_face_values.JxW(q);
						scratch.u2_v2(i,j) += -0.5*(-phi_grad_i_n*phi_j_n-phi_grad_j_n*phi_i_n)*normal_v*
								scratch.fe_face_values.JxW(q)+
								penalty*(phi_i_n*phi_j_n)*scratch.fe_face_values.JxW(q);
					}
				}
			}//endforq
			copy.neighbor_u_v[face_n][0] = scratch.u1_v1;
			copy.neighbor_u_v[face_n][1] = scratch.u2_v2;
			copy.neighbor_u_v[face_n][2] = scratch.u1_v2;
			copy.neighbor_u_v[face_n][3] = scratch.u2_v1;

			neighbor->get_dof_indices(copy.neighbor_dof_indices[face_n]);
		}//endif
	}//endface
		}
template <int dim>
void Assembler_s<dim>::copy_cell_diff(const Copy_diff &copy)
{
	for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
			diffusion_matrix_s.add(copy.local_dof_indices[i],
					copy.local_dof_indices[j],
					copy.local_matrix(i,j));
	for (unsigned int face_n=0; face_n<4; ++face_n)
	{
		if (copy.sides[face_n])
			for (unsigned int i=0; i<dofs_per_cell; ++i)
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
					diffusion_matrix_s.add(copy.local_dof_indices[i],
							copy.local_dof_indices[j],
							0.5*copy.neighbor_u_v[face_n][0](i,j));
					diffusion_matrix_s.add(copy.neighbor_dof_indices[face_n][i],
							copy.neighbor_dof_indices[face_n][j],
							0.5*copy.neighbor_u_v[face_n][1](i,j));
					diffusion_matrix_s.add(copy.neighbor_dof_indices[face_n][i],
							copy.local_dof_indices[j],
							0.5*copy.neighbor_u_v[face_n][2](i,j));
					diffusion_matrix_s.add(copy.local_dof_indices[i],
							copy.neighbor_dof_indices[face_n][j],
							0.5*copy.neighbor_u_v[face_n][3](i,j));
				}
	}
}
template <int dim>
void Assembler_s<dim>::assemble_diffusion_matrix_s()
{
	QGauss<dim> quadrature_formula(fe_degree+1);
	QGauss<dim-1> face_quadrature_formula(fe_degree+1);
	WorkStream::run(dof_handler->begin_active(),
			dof_handler->end(),
			std::bind(&Assembler_s<dim>::assemble_cell_diff,
					this,
					std::placeholders::_1,
					std::placeholders::_2,
					std::placeholders::_3),
					std::bind(&Assembler_s<dim>::copy_cell_diff,
							this,
							std::placeholders::_1),
							Scratch_diff(*fe,quadrature_formula,
									face_quadrature_formula),
									Copy_diff());
}
template <int dim>
Assembler_s<dim>::Scratch_source::Scratch_source(const FiniteElement<dim> &fe,
		const Quadrature<dim> &quadrature)
:
		fe_values(fe, quadrature,
				update_values | update_gradients|
				update_quadrature_points| update_JxW_values),
				function_values(4)
{
	for (unsigned int i=0; i<4; ++i)
		function_values[i].resize(quadrature.size());
}

template <int dim>
Assembler_s<dim>::Scratch_source::Scratch_source(const Scratch_source &scratch)
:
fe_values(scratch.fe_values.get_fe(),
		scratch.fe_values.get_quadrature(),
		scratch.fe_values.get_update_flags()),
		function_values(scratch.function_values)
{}

template <int dim>
void Assembler_s<dim>::assemble_cell_source(const typename DoFHandler<dim>::active_cell_iterator &cell,
		Scratch_source &scratch, Copy_source &copy,
		const std::vector<Vector<double> > &functions)
		{
	const unsigned int n_q_points =scratch.fe_values.get_quadrature().size();
	copy.local_dof_indices.resize(dofs_per_cell);
	cell->get_dof_indices(copy.local_dof_indices);
	copy.source_values.resize(2);

	scratch.fe_values.reinit(cell);
	for (unsigned int i=0; i<2; ++i)
	{
		scratch.fe_values.get_function_values(functions[i],
				scratch.function_values[i]);
		copy.source_values[i].reinit(dofs_per_cell);
	}

	for (unsigned int q=0; q<n_q_points; ++q)
	{
		const double M0_source = biology_s->M0_source(scratch.function_values[0][q],
				scratch.function_values[1][q]);
		const double S_source = biology_s->S_source(scratch.function_values[0][q],
				scratch.function_values[1][q]);

		for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
			copy.source_values[0][i] += M0_source*scratch.fe_values.shape_value(i,q)*
					scratch.fe_values.JxW(q);
			copy.source_values[1][i] += S_source*scratch.fe_values.shape_value(i,q)*
					scratch.fe_values.JxW(q);

		}
	}
		}

template <int dim>
void Assembler_s<dim>::copy_cell_source(const Copy_source &copy,
		std::vector<Vector<double> > &source)
		{
	for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<2; ++j)////////////////////////////////////////////////////////////////////////////////////// j 3?
			source[j][copy.local_dof_indices[i]] = copy.source_values[j][i];
		}

template <int dim>
void Assembler_s<dim>::assemble_source(const std::vector<Vector<double> > &functions,
		std::vector<Vector<double> > &source)
		{
	QGauss<dim> quadrature_formula(fe_degree+1);
	WorkStream::run(dof_handler->begin_active(),
			dof_handler->end(),
			std::bind(&Assembler_s<dim>::assemble_cell_source,
					this,
					std::placeholders::_1,
					std::placeholders::_2,
					std::placeholders::_3,
					std::ref(functions)),
					std::bind(&Assembler_s<dim>::copy_cell_source,
							this,
							std::placeholders::_1,
							std::ref(source)),
							Scratch_source(*fe,quadrature_formula),
							Copy_source());
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <int dim>
void Assembler_s<dim>::compute_p_d(const std::vector<double> &points,
		const std::vector<double> &func_values,
		const std::vector<bool> &sides,
		std::vector<double>  &derivatives) const
		{


	if (sides[0]&&sides[1])
	{
		derivatives[2] = (func_values[1]-func_values[0])/(points[1]-points[0]);
		derivatives[0] = (func_values[2]-func_values[0])/(points[2]-points[0]);
		derivatives[1] = (func_values[1]-func_values[2])/(points[1]-points[2]);
	}
	else
	{

		if(sides[0])
		{
			derivatives[2] = (func_values[2]-func_values[0])/(points[2]-points[0]);
			derivatives[0] = (func_values[2]-func_values[0])/(points[2]-points[0]);
		}
		if(sides[1])
		{
			derivatives[2] = (func_values[1]-func_values[2])/(points[1]-points[2]);
			derivatives[1] = (func_values[1]-func_values[2])/(points[1]-points[2]);
		}
	}

		}
template <int dim>
void Assembler_s<dim>::compute_drift_field(const std::vector<Point<dim> > &points,
		const std::vector<Tensor<2,dim> > &DT_values,
		const std::vector<bool> &sides,
		std::vector<Tensor<1,dim> > &drifts) const
		{

	std::vector<double> x_points(3);
	x_points[0] = points[0][0];
	x_points[1] = points[1][0];
	x_points[2] = points[4][0];
	std::vector<bool> sides_x(2);
	sides_x[0] = sides[0];
	sides_x[1] = sides[1];

	std::vector<double> y_points(3);
	y_points[0] = points[2][1];
	y_points[1] = points[3][1];
	y_points[2] = points[4][1];
	std::vector<bool> sides_y(2);
	sides_y[0] = sides[2];
	sides_y[1] = sides[3];

	Tensor<2,dim> dx_DT;
	Tensor<2,dim> dy_DT;
	Tensor<2,dim> dx_DT_L;
	Tensor<2,dim> dx_DT_R;
	Tensor<2,dim> dy_DT_D;
	Tensor<2,dim> dy_DT_U;

	for (unsigned int i=0; i<dim; ++i)
	{
		for(unsigned int j=0; j<dim; ++j)
		{

			std::vector<double> DTij_x(3);
			DTij_x[0] = DT_values[0][i][j];
			DTij_x[1] = DT_values[1][i][j];
			DTij_x[2] = DT_values[4][i][j];

			std::vector<double> derivative_x(3);
			Assembler_s::compute_p_d(x_points,DTij_x,sides_x,derivative_x);

			std::vector<double> DTij_y(3);
			DTij_y[0] = DT_values[2][i][j];
			DTij_y[1] = DT_values[3][i][j];
			DTij_y[2] = DT_values[4][i][j];

			std::vector<double> derivative_y(3);
			Assembler_s::compute_p_d(y_points,DTij_y,sides_y,derivative_y);

			dx_DT[i][j] = derivative_x[2];
			dy_DT[i][j] = derivative_y[2];

			dx_DT_L[i][j] = derivative_x[0];
			dx_DT_R[i][j] = derivative_x[1];

			dy_DT_D[i][j] = derivative_y[0];
			dy_DT_U[i][j] = derivative_y[1];
		}
	}


	Tensor<1,dim> div_DT;
	div_DT[0] = dx_DT[0][0]+dy_DT[0][1];
	div_DT[1] = dx_DT[1][0]+dy_DT[1][1];

	drifts[4] = div_DT;

	Tensor<1,dim> div_DT_L;
	div_DT_L[0] = dx_DT_L[0][0]+dy_DT[0][1];
	div_DT_L[1] = dx_DT_L[1][0]+dy_DT[1][1];

	drifts[0] = div_DT_L;

	Tensor<1,dim> div_DT_R;
	div_DT_R[0] = dx_DT_R[0][0]+dy_DT[0][1];
	div_DT_R[1] = dx_DT_R[1][0]+dy_DT[1][1];

	drifts[1] = div_DT_R;

	Tensor<1,dim> div_DT_D;
	div_DT_D[0] = dx_DT[0][0]+dy_DT_D[0][1];
	div_DT_D[1] = dx_DT[1][0]+dy_DT_D[1][1];

	drifts[2] = div_DT_D;

	Tensor<1,dim> div_DT_U;
	div_DT_U[0] = dx_DT[0][0]+dy_DT_U[0][1];
	div_DT_U[1] = dx_DT[1][0]+dy_DT_U[1][1];

	drifts[3] = div_DT_U;
		}

template <int dim>
void Assembler_s<dim>::compute_sub_func(const double &S, double &g)
{
	//const double g1 = lambda_1/(k_p*S+k_m+lambda_0);
	//const double g2 = k_p*R_0/(k_p*S+k_m) - pow(k_p,2)*S*R_0/pow(k_p*S+k_m,2);

	const double g1 = lambda_1/(S+k_new+lambda_0);
	const double g2 = k_new/pow(S+k_new,2);

	g = g1*g2;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <int dim>
Assembler_s<dim>::Scratch_diff2::Scratch_diff2(const FiniteElement<dim> &fe,
		const Quadrature<dim> &quadrature,
		const Quadrature<dim-1> &face_quadrature)
:
		fe_values(fe, quadrature,
				update_values | update_gradients|
				update_quadrature_points| update_JxW_values),
				fe_face_values(fe, face_quadrature,
						update_values| update_gradients| update_normal_vectors|
						update_quadrature_points | update_JxW_values),
						fe_face_values_neighbor(fe, face_quadrature,
								update_values| update_gradients| update_normal_vectors|
								update_quadrature_points | update_JxW_values),
								u1_v1(fe.dofs_per_cell, fe.dofs_per_cell),
								u2_v2(fe.dofs_per_cell, fe.dofs_per_cell),
								u1_v2(fe.dofs_per_cell, fe.dofs_per_cell),
								u2_v1(fe.dofs_per_cell, fe.dofs_per_cell)
{}

template <int dim>
Assembler_s<dim>::Scratch_diff2::Scratch_diff2(const Scratch_diff2 &scratch)
:
fe_values(scratch.fe_values.get_fe(),
		scratch.fe_values.get_quadrature(),
		scratch.fe_values.get_update_flags()),
		fe_face_values(scratch.fe_face_values.get_fe(),
				scratch.fe_face_values.get_quadrature(),
				scratch.fe_face_values.get_update_flags()),
				fe_face_values_neighbor(scratch.fe_face_values_neighbor.get_fe(),
						scratch.fe_face_values_neighbor.get_quadrature(),
						scratch.fe_face_values_neighbor.get_update_flags()),
						u1_v1(scratch.u1_v1),
						u2_v2(scratch.u2_v2),
						u1_v2(scratch.u1_v2),
						u2_v1(scratch.u2_v1)
{}

template <int dim>
void Assembler_s<dim>::assemble_cell_diff2(const typename DoFHandler<dim>::active_cell_iterator &cell,
		Scratch_diff2 &scratch, Copy_diff2 &copy)
		{
	const unsigned int n_q_points = scratch.fe_values.get_quadrature().size();
	const unsigned int n_face_q_points = scratch.fe_face_values.get_quadrature().size();

	copy.local_matrix.reinit(dofs_per_cell, dofs_per_cell);
	copy.neighbor_u_v.resize(4);
	copy.neighbor_dof_indices.resize(4);
	copy.sides.resize(4);
	copy.local_dof_indices.resize(dofs_per_cell);
	cell->get_dof_indices(copy.local_dof_indices);
	//const double edge_length = sqrt(0.5)*cell->diameter();
	/////////////////////////////////////////////////////////
	Tensor<2,dim> DT_in;
	Tensor<2,dim> DT_out;
	Point<2> pt;
	pt[0]=100;
	pt[1]=150;

	DT_in = diffusionTensor.value(pt);

	std::vector<Tensor<2,dim> > DT_values(5);
	std::vector<Tensor<1,dim> > drifts(5);
	std::vector<bool> sides(4);
	std::vector<Point<dim> > cell_Points(5);
	cell_Points[4] = cell->center();
	DT_values[4] = DT_in;
	for (unsigned int face_n=0; face_n<GeometryInfo<dim>::faces_per_cell;
			++face_n)
	{
		if (!cell->at_boundary(face_n))
		{
			sides[face_n] = true;
			typename DoFHandler<dim>::active_cell_iterator
			neighbor = cell->neighbor(face_n);
			cell_Points[face_n] = neighbor->center();
			DT_values[face_n] = DT_in;
		}
		else
		{
			sides[face_n] = false;
		}
	}
	DT_in *= T_scale;
	//	std::cout << " DTI "<< DT_in<<std::endl;
	//   std::cout<< "..................."<<std::endl;
	//    std::cout<< T_scale<<std::endl;

	for (unsigned int i=0; i<5; ++i)
		drifts[i] *= 2*T_scale;

	//	std::cout<< "drifts "<<"  " <<drifts[0] << "  " <<drifts[1]<<"  " <<drifts[2] << "  " <<
	//			drifts[3]<<"  " <<drifts[4]<<std::endl;
	//	std::cout<< "..................."<<std::endl;
	////////////////////////////////////////////////////////////////
	scratch.fe_values.reinit(cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
			for (unsigned int q=0; q<n_q_points; ++q)
				copy.local_matrix(i,j) += DT_in*scratch.fe_values.shape_grad(i,q)*/////////////////////
				scratch.fe_values.shape_grad(j,q)*scratch.fe_values.JxW(q);
	////////////////////////////////////////////////////////////////
	for (unsigned int face_n=0; face_n<4; ++face_n)
	{
		scratch.u1_v1 = 0;
		scratch.u1_v2 = 0;
		scratch.u2_v1 = 0;
		scratch.u2_v2 = 0;

		copy.neighbor_dof_indices[face_n].resize(dofs_per_cell);
		copy.neighbor_u_v[face_n].resize(4);
		for (unsigned int i=0; i<4; ++i)
			copy.neighbor_u_v[face_n][i].reinit(dofs_per_cell, dofs_per_cell);
		copy.sides[face_n] = false;

		if (!cell->at_boundary(face_n))
		{
			DT_out = DT_values[face_n];
			DT_out *= T_scale;


			Tensor<2,dim> DT_int = DT_out + DT_in;
			DT_int *= 0.5;

			scratch.fe_face_values.reinit(cell,face_n);
			typename DoFHandler<dim>::active_cell_iterator
			neighbor = cell->neighbor(face_n);
			unsigned int neighbor_face = cell->neighbor_of_neighbor(face_n);
			scratch.fe_face_values_neighbor.reinit(neighbor,neighbor_face);
			copy.sides[face_n] = true;
			const double penalty = 100;//PENALTY HERE
			for (unsigned int q=0; q<n_face_q_points; ++q)
			{
				const Tensor<1,dim> normal_v = scratch.fe_face_values.normal_vector(q);
				for (unsigned int i=0; i<dofs_per_cell; ++i)
				{
					const Tensor<1,dim> phi_grad_i = scratch.fe_face_values.shape_grad(i,q);
					const double phi_i = scratch.fe_face_values.shape_value(i,q);
					const Tensor<1,dim> phi_grad_i_n = scratch.fe_face_values_neighbor.shape_grad(i,q);
					const double phi_i_n = scratch.fe_face_values_neighbor.shape_value(i,q);
					for (unsigned int j=0; j<dofs_per_cell; ++j)
					{
						const Tensor<1,dim> phi_grad_j = scratch.fe_face_values.shape_grad(j,q);
						const double phi_j = scratch.fe_face_values.shape_value(j,q);
						const Tensor<1,dim> phi_grad_j_n = scratch.fe_face_values_neighbor.shape_grad(j,q);
						const double phi_j_n = scratch.fe_face_values_neighbor.shape_value(j,q);

						scratch.u1_v1(i,j) +=(-0.5*((DT_out*(phi_i*phi_grad_j+phi_j*phi_grad_i))*normal_v)+
								penalty*phi_i*phi_j)*scratch.fe_face_values.JxW(q);


						scratch.u1_v2(i,j) += (0.5*(DT_in*phi_grad_j)*normal_v*phi_i_n-0.5*(DT_out*phi_grad_i_n)*normal_v*phi_j-
								penalty*phi_i_n*phi_j)*scratch.fe_face_values.JxW(q);


						scratch.u2_v1(i,j) += (-0.5*(DT_out*phi_grad_j_n)*normal_v*phi_i+0.5*(DT_in*phi_grad_i)*normal_v*phi_j_n-
								penalty*phi_i*phi_j_n)*scratch.fe_face_values.JxW(q);


						scratch.u2_v2(i,j) += (0.5*((DT_out*(phi_i_n*phi_grad_j_n+phi_j_n*phi_grad_i_n))*normal_v)+
								penalty*phi_i_n*phi_j_n)*scratch.fe_face_values.JxW(q);

					}
				}
			}//endforq
			copy.neighbor_u_v[face_n][0] = scratch.u1_v1;
			copy.neighbor_u_v[face_n][1] = scratch.u2_v2;
			copy.neighbor_u_v[face_n][2] = scratch.u1_v2;
			copy.neighbor_u_v[face_n][3] = scratch.u2_v1;

			neighbor->get_dof_indices(copy.neighbor_dof_indices[face_n]);
		}//endif
	}//endface

		}
template <int dim>
void Assembler_s<dim>::copy_cell_diff2(const Copy_diff2 &copy)
{
	for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
			diffusion_matrix_m.add(copy.local_dof_indices[i],
					copy.local_dof_indices[j],
					copy.local_matrix(i,j));
	for (unsigned int face_n=0; face_n<4; ++face_n)
	{
		if (copy.sides[face_n])
			for (unsigned int i=0; i<dofs_per_cell; ++i)
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
					diffusion_matrix_m.add(copy.local_dof_indices[i],
							copy.local_dof_indices[j],
							0.5*copy.neighbor_u_v[face_n][0](i,j));
					diffusion_matrix_m.add(copy.neighbor_dof_indices[face_n][i],
							copy.neighbor_dof_indices[face_n][j],
							0.5*copy.neighbor_u_v[face_n][1](i,j));
					diffusion_matrix_m.add(copy.neighbor_dof_indices[face_n][i],
							copy.local_dof_indices[j],
							0.5*copy.neighbor_u_v[face_n][2](i,j));
					diffusion_matrix_m.add(copy.local_dof_indices[i],
							copy.neighbor_dof_indices[face_n][j],
							0.5*copy.neighbor_u_v[face_n][3](i,j));
				}
	}
}
template <int dim>
void Assembler_s<dim>::assemble_diffusion_matrix_m()
{
	QGauss<dim> quadrature_formula(fe_degree+1);
	QGauss<dim-1> face_quadrature_formula(fe_degree+1);
	WorkStream::run(dof_handler->begin_active(),
			dof_handler->end(),
			std::bind(&Assembler_s<dim>::assemble_cell_diff2,
					this,
					std::placeholders::_1,
					std::placeholders::_2,
					std::placeholders::_3),
					std::bind(&Assembler_s<dim>::copy_cell_diff2,
							this,
							std::placeholders::_1),
							Scratch_diff2(*fe,quadrature_formula,
									face_quadrature_formula),
									Copy_diff2());
}
//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

template <int dim>
Assembler_s<dim>::Scratch_adv::Scratch_adv(const FiniteElement<dim> &fe,
		const Quadrature<dim> &quadrature,
		const Quadrature<dim-1> &face_quadrature)
:
		fe_values(fe, quadrature,
				update_values | update_gradients|
				update_quadrature_points| update_JxW_values),
				fe_values_neighbor(fe, quadrature,
						update_values | update_gradients|
						update_quadrature_points| update_JxW_values),
						fe_face_values(fe, face_quadrature,
								update_values| update_gradients| update_normal_vectors|
								update_quadrature_points | update_JxW_values),
								fe_face_values_neighbor(fe, face_quadrature,
										update_values| update_gradients| update_normal_vectors|
										update_quadrature_points | update_JxW_values),
										u1_v1(fe.dofs_per_cell, fe.dofs_per_cell),
										u2_v2(fe.dofs_per_cell, fe.dofs_per_cell),
										u1_v2(fe.dofs_per_cell, fe.dofs_per_cell),
										u2_v1(fe.dofs_per_cell, fe.dofs_per_cell)
{}

template <int dim>
Assembler_s<dim>::Scratch_adv::Scratch_adv(const Scratch_adv &scratch)
:
fe_values(scratch.fe_values.get_fe(),
		scratch.fe_values.get_quadrature(),
		scratch.fe_values.get_update_flags()),
		fe_values_neighbor(scratch.fe_values_neighbor.get_fe(),
				scratch.fe_values_neighbor.get_quadrature(),
				scratch.fe_values_neighbor.get_update_flags()),
				fe_face_values(scratch.fe_face_values.get_fe(),
						scratch.fe_face_values.get_quadrature(),
						scratch.fe_face_values.get_update_flags()),
						fe_face_values_neighbor(scratch.fe_face_values_neighbor.get_fe(),
								scratch.fe_face_values_neighbor.get_quadrature(),
								scratch.fe_face_values_neighbor.get_update_flags()),
								u1_v1(scratch.u1_v1),
								u2_v2(scratch.u2_v2),
								u1_v2(scratch.u1_v2),
								u2_v1(scratch.u2_v1)
{}

template <int dim>
void Assembler_s<dim>::assemble_cell_adv(const typename DoFHandler<dim>::active_cell_iterator &cell,
		Scratch_adv &scratch, Copy_adv &copy,
		const Vector<double> &h)

		{
	const unsigned int n_q_points = scratch.fe_values.get_quadrature().size();
	const unsigned int n_face_q_points = scratch.fe_face_values.get_quadrature().size();
	const double volume = cell->measure();
	const double delta = 1./std::sqrt(2)*cell->diameter();
	copy.local_matrix.reinit(dofs_per_cell, dofs_per_cell);
	copy.neighbor_u_v.resize(4);
	copy.neighbor_dof_indices.resize(4);
	copy.sides.resize(4);
	copy.local_dof_indices.resize(dofs_per_cell);
	cell->get_dof_indices(copy.local_dof_indices);

	scratch.fe_values.reinit(cell);
	std::vector<double> h_values(4);
	std::vector<double> h_q_values(n_q_points);
	std::vector<Tensor<1,dim> > grad_h(5);
	for (unsigned int face_n=0; face_n<4; ++face_n)
	{
		copy.sides[face_n] = false;
		if (!cell->at_boundary(face_n))
		{
			copy.sides[face_n] = true;
			typename DoFHandler<dim>::active_cell_iterator
			neighbor = cell->neighbor(face_n);
			scratch.fe_values_neighbor.reinit(neighbor);
			scratch.fe_values_neighbor.get_function_values(h,h_q_values);

			for (unsigned int q=0; q<n_q_points; ++q)
			{
				h_values[face_n] +=1./volume*h_q_values[q]*scratch.fe_values_neighbor.JxW(q);

			}
		}
	}


	double cell_h=0.0;
	std::vector<double> cell_h_values(n_q_points);
	scratch.fe_values.get_function_values(h,cell_h_values);

	for (unsigned int q=0; q<n_q_points; ++q)
	{
		cell_h += 1./volume*cell_h_values[q]*scratch.fe_values.JxW(q);


	}

	if (copy.sides[0])
	{
		grad_h[0][0]= (cell_h-h_values[0])/delta;

	}
	if (copy.sides[1])
	{
		grad_h[1][0] = (h_values[1]-cell_h)/delta;

	}
	if (copy.sides[2])
	{
		grad_h[2][1] = (cell_h-h_values[2])/delta;

	}
	if (copy.sides[3])
	{
		grad_h[3][1] = (h_values[3]-cell_h)/delta;

	}

	if (copy.sides[0]&&copy.sides[1])
	{
		grad_h[4][0] = 0.5*(grad_h[0][0]+grad_h[1][0]);

	}
	else if (copy.sides[0])
	{
		grad_h[4][0] = grad_h[0][0];

	}
	else
	{
		grad_h[4][0] = grad_h[1][0];

	}

	if (copy.sides[2]&&copy.sides[3])
	{
		grad_h[4][1] = 0.5*(grad_h[2][1]+grad_h[3][1]);

	}
	else if (copy.sides[2])
	{
		grad_h[4][1] = grad_h[2][1];

	}
	else
	{
		grad_h[4][1] = grad_h[3][1];

	}

	Tensor<2,dim> DT_in;
	Tensor<2,dim> DT_out;
	Point<2> pt;
	pt[0]=100;
	pt[1]=150;
	DT_in = diffusionTensor.value(pt);
	std::vector<Tensor<2,dim> > DT_values(5);
	std::vector<Tensor<1,dim> > drifts(5);
	DT_values[4] = DT_in;
	std::vector<bool> sides(4);
	std::vector<Point<dim> > cell_Points(5);
	cell_Points[4] = cell->center();
	double g_func_n;

	for (unsigned int face_n=0; face_n<GeometryInfo<dim>::faces_per_cell;
			++face_n)
	{
		if (!cell->at_boundary(face_n))
		{
			sides[face_n] = true;
			typename DoFHandler<dim>::active_cell_iterator
			neighbor = cell->neighbor(face_n);
			cell_Points[face_n] = neighbor->center();
			DT_values[face_n] = DT_in;

		}
		else
		{
			sides[face_n] = false;
		}
	}


	DT_in *= T_scale;

	Assembler_s::compute_drift_field(cell_Points,DT_values,sides,drifts);

	for (unsigned int i=0; i<5; ++i)
		drifts[i] *= 2*T_scale;

	//std::cout<< "drifts "<<"  " <<drifts[0] << "  " <<drifts[1]<<"  " <<drifts[2] << "  " <<
	//		drifts[3]<<"  " <<drifts[4]<<std::endl;
	//std::cout<< "..................."<<std::endl;
	//std::cout<<"DT_in  : "<<DT_in<<std::endl;


	Assembler_s::compute_sub_func(cell_h, g_func_n);
	Tensor<1,dim> velocity = (g_func_n*DT_in*grad_h[4]+drifts[4]);
	//Tensor<1,dim> velocity;
	//velocity[0] = 0;
	//velocity[1] = 0;


//	    std::cout << " velocity "<< velocity<<std::endl;
//	    std::cout<< "..................."<<std::endl;
	for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
			for (unsigned int q=0; q<n_q_points; ++q)
			{
				{
					{

						copy.local_matrix(i,j) += -velocity*scratch.fe_values.shape_grad(i,q)*
								scratch.fe_values.shape_value(j,q)*scratch.fe_values.JxW(q);
					}
				}
			}
	/////////////////////////////////////


	for (unsigned int face_n=0; face_n<4; ++face_n)
	{
		scratch.u1_v1 = 0;
		scratch.u1_v2 = 0;
		scratch.u2_v1 = 0;
		scratch.u2_v2 = 0;

		copy.neighbor_dof_indices[face_n].resize(dofs_per_cell);
		copy.neighbor_u_v[face_n].resize(4);
		for (unsigned int i=0; i<4; ++i)
			copy.neighbor_u_v[face_n][i].reinit(dofs_per_cell, dofs_per_cell);

		if (!cell->at_boundary(face_n))
		{
			DT_out = DT_values[face_n];
			DT_out *= T_scale;
			Tensor<2,dim> DT_int = DT_out + DT_in;
			DT_int *= 0.5;

			double g_func_n2;
			Assembler_s::compute_sub_func(h_values[face_n], g_func_n2);
			scratch.fe_face_values.reinit(cell,face_n);
			typename DoFHandler<dim>::active_cell_iterator
			neighbor = cell->neighbor(face_n);
			unsigned int neighbor_face = cell->neighbor_of_neighbor(face_n);
			scratch.fe_face_values_neighbor.reinit(neighbor,neighbor_face);



			for (unsigned int q=0; q<n_face_q_points; ++q)
			{

				const Tensor<1,dim> normal = scratch.fe_face_values.normal_vector(q); //////////////////////////////////////////



				for (unsigned int i=0; i<dofs_per_cell; ++i)
					for (unsigned int j=0; j<dofs_per_cell; ++j)
					{

						const double flux = -normal*(g_func_n2*DT_int*grad_h[face_n]+
							drifts[face_n]);
						//const double flux = -normal*velocity;

						if (flux>=0)
						{

							scratch.u1_v1(i,j) += flux*scratch.fe_face_values.shape_value(i,q)*
									scratch.fe_face_values.shape_value(j,q)*scratch.fe_face_values.JxW(q);
							scratch.u1_v2(i,j) += -flux*scratch.fe_face_values_neighbor.shape_value(i,q)*
									scratch.fe_face_values.shape_value(j,q)*scratch.fe_face_values.JxW(q);
						}

						else
						{
							scratch.u2_v1(i,j) += flux*scratch.fe_face_values.shape_value(i,q)*
									scratch.fe_face_values_neighbor.shape_value(j,q)*
									scratch.fe_face_values.JxW(q);
							scratch.u2_v2(i,j) += -flux*scratch.fe_face_values_neighbor.shape_value(i,q)*
									scratch.fe_face_values_neighbor.shape_value(j,q)*
									scratch.fe_face_values.JxW(q);
						}
					}
			}
			copy.neighbor_u_v[face_n][0] = scratch.u1_v1;
			copy.neighbor_u_v[face_n][1] = scratch.u2_v2;
			copy.neighbor_u_v[face_n][2] = scratch.u1_v2;
			copy.neighbor_u_v[face_n][3] = scratch.u2_v1;

			neighbor->get_dof_indices(copy.neighbor_dof_indices[face_n]);
		}//endif
	}//endface
		}

template <int dim>
void Assembler_s<dim>::copy_cell_adv(const Copy_adv &copy)
{
	for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
			advection_matrix.add(copy.local_dof_indices[i],
					copy.local_dof_indices[j],
					copy.local_matrix(i,j));
	for (unsigned int face_n=0; face_n<4; ++face_n)
	{
		if (copy.sides[face_n])
			for (unsigned int i=0; i<dofs_per_cell; ++i)
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
					advection_matrix.add(copy.local_dof_indices[i],
							copy.local_dof_indices[j],
							0.5*copy.neighbor_u_v[face_n][0](i,j));
					advection_matrix.add(copy.neighbor_dof_indices[face_n][i],
							copy.neighbor_dof_indices[face_n][j],
							0.5*copy.neighbor_u_v[face_n][1](i,j));
					advection_matrix.add(copy.neighbor_dof_indices[face_n][i],
							copy.local_dof_indices[j],
							0.5*copy.neighbor_u_v[face_n][2](i,j));
					advection_matrix.add(copy.local_dof_indices[i],
							copy.neighbor_dof_indices[face_n][j],
							0.5*copy.neighbor_u_v[face_n][3](i,j));
				}
	}
}

template <int dim>
void Assembler_s<dim>::assemble_advection_matrix(const Vector<double> &h)

{
	advection_matrix.reinit(flux_sparsity);
	QGauss<dim> quadrature_formula(fe_degree+1);
	QGauss<dim-1> face_quadrature_formula(fe_degree+1);

	WorkStream::run(dof_handler->begin_active(),
			dof_handler->end(),
			std::bind(&Assembler_s<dim>::assemble_cell_adv,
					this,
					std::placeholders::_1,
					std::placeholders::_2,
					std::placeholders::_3,
					std::ref(h)),

					std::bind(&Assembler_s<dim>::copy_cell_adv,
							this,
							std::placeholders::_1),
							Scratch_adv(*fe,quadrature_formula,
									face_quadrature_formula),
									Copy_adv());
}

template <int dim>
void Assembler_s<dim>::initial_values(std::vector<Vector<double> > &initial) const
{
	QGauss<dim> quadrature_formula(fe_degree+1);
	FEValues<dim> fe_values(*fe, quadrature_formula,
			update_values| update_quadrature_points|
			update_JxW_values);
	const unsigned int n_q_points = quadrature_formula.size();
	std::vector<Point<dim> > quadrature_points(n_q_points);
	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	std::vector<Vector<double> > func_rhs(2);
	for (unsigned int i=0; i<2; ++i)
		func_rhs[i].reinit(n_dofs);
	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler->begin_active(),
	endc = dof_handler->end();
	for (; cell!=endc; ++cell)
	{
		fe_values.reinit(cell);
		cell->get_dof_indices(local_dof_indices);
		quadrature_points = fe_values.get_quadrature_points();
		for (unsigned int q=0; q<n_q_points; ++q)
		{
			const double M0 = biology_s->initial_M0(quadrature_points[q]);
			const double S = biology_s->initial_S(quadrature_points[q]);


			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				func_rhs[0](local_dof_indices[i]) += M0*fe_values.shape_value(i,q)*
						fe_values.JxW(q);
				func_rhs[1](local_dof_indices[i]) += S*fe_values.shape_value(i,q)*
						fe_values.JxW(q);

			}
		}
	}
	for (unsigned int i=0; i<2; ++i)
		solve_mass(func_rhs[i], initial[i]);
}

template <int dim>
void Assembler_s<dim>::solve_mass(const Vector<double> &rhs,
		Vector<double> &solution) const
		{
	SolverControl solver_control(5000,1e-12);
	SolverCG<> solver(solver_control);
	solver.solve(mass_mass_matrix, solution, rhs, mass_precond);
		}

#endif
