% this script calculate the eigen values of J3 for a given c

clear all;
clc;
close all;

%% dimensionalized data
M_max = 0.8;
S_max = 10^(-6.4);
alpha_old = 1e-11;
beta_old = 1e-9; 
%% non-dimensionalized data
alpha = (alpha_old*S_max)/(beta_old*M_max);
D = 1/100;
mu_old = 0.2/(60*60*24);
mu0 = (mu_old*S_max)/(beta_old*M_max);
%%
c_old3 = 4*sqrt(D*mu0);
% c_old3 = 1e+10;
g3 = function_g(1);
a1 = c_old3 + c_old3/D;
a2 = ((c_old3^2)/D) - alpha*g3 - alpha;
a3 = -c_old3*alpha/D;
a4 = alpha*(mu0/D)*(1-alpha);

%eigen values
p = [1 a1 a2 a3 a4];
eigJ2_old = roots(p);
