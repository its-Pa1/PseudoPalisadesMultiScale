//SGN
#ifndef BIOLOGY_HYPO_MICRO_HH
#define BIOLOGY_HYPO_MICRO_HH

#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <numeric>
using namespace dealii;

class Biology_s
{
public:
	Biology_s();
	//double M0_source(const double m0,const double s,const double &Qx) const;
	double M0_source(const double m0,const double s) const;
	double S_source(const double m0, const double s) const;

	double initial_M0(const Point<2> &p) const;
	double initial_S(const Point<2> &p) const;

	const double alpha = 0.001;
	const double D_s = 100;
	const double D_m = 1; //
	const double beta = 1;
private:

	const unsigned int dim = 2;
	double cell_size;

};

Biology_s::Biology_s()
{}


//double Biology_s::M0_source(const double m0,const double s,const double &Qx) const/////////////////////////// here i was
double Biology_s::M0_source(const double m0,const double s) const
{
	//return (0.01-s)*Qx*m0;//third
	return ((1e-7)-s)*(1-m0)*m0;//first
	//return 0.1*m0;
	//return ;
}

double Biology_s::S_source(const double m0, const double s) const
{
	return beta*m0-alpha*s;////////////////////////////////////////////////////////////
}

double Biology_s::initial_M0(const Point<2> &p) const
{

	const double radius = 50; const double radius_2 = 40; const double radius_3 = 30;
	const double sigma = 20;   const double sigma_2 = 15;   const double sigma_3 = 10;

	const double centerx = 500;  const double centerx_2 = 500;  const double centerx_3 = 600;
	const double centery = 500;  const double centery_2 = 600;  const double centery_3 = 500;

	const double dx =  pow(p[0]-centerx,2);
	const double dy =  pow(p[1]-centery,2);
	const double distance = sqrt(dx+dy );

	const double dx_2 =  pow(p[0]-centerx_2,2);
	const double dy_2 =  pow(p[1]-centery_2,2);
	const double distance_2 = sqrt(dx_2+dy_2 );

	const double dx_3 =  pow(p[0]-centerx_3,2);
	const double dy_3 =  pow(p[1]-centery_3,2);
	const double distance_3 = sqrt(dx_3+dy_3 );

	double f_val = 0.5*exp(-(dx+dy)/(2*pow(sigma,2)));
	double f_val_2 = 0.5*exp(-(dx_2+dy_2)/(2*pow(sigma_2,2)));
	double f_val_3 = 0.5*exp(-(dx_3+dy_3)/(2*pow(sigma_3,2)));
	if (distance<radius)
		return f_val;
		else if (distance_2<radius_2)
			return f_val_2;
		else if (distance_3<radius_3)
			return f_val_3;
	else
		return 0;
}

double Biology_s::initial_S(const Point<2> &p) const
{
	const double radius = 30;     const double radius_2 = 20;  const double radius_3 = 15;
	const double sigma = 10;       const double sigma_2 = 8;    const double sigma_3 = 7;

	const double centerx = 500;  const double centerx_2 = 500;  const double centerx_3 = 600;
	const double centery = 500;  const double centery_2 = 600;  const double centery_3 = 500;

	const double dx =  pow(p[0]-centerx,2);
	const double dy =  pow(p[1]-centery,2);
	const double distance = sqrt(dx+dy );

	const double dx_2 =  pow(p[0]-centerx_2,2);
	const double dy_2 =  pow(p[1]-centery_2,2);
	const double distance_2 = sqrt(dx_2+dy_2 );

	const double dx_3 =  pow(p[0]-centerx_3,2);
	const double dy_3 =  pow(p[1]-centery_3,2);
	const double distance_3 = sqrt(dx_3+dy_3);

	double f_val = (1e-6)*exp(-(dx+dy)/(2*pow(sigma,2)));
	double f_val_2 = (1e-6)*exp(-(dx_2+dy_2)/(2*pow(sigma_2,2)));
	double f_val_3 = (1e-6)*exp(-(dx_3+dy_3)/(2*pow(sigma_3,2)));
	if (distance<radius)
		return f_val;
		else if (distance_2<radius_2)
			return f_val_2;
		else if (distance_3<radius_3)
			return f_val_3;
	else
		return 0;
}


#endif
