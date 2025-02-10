% This script computes the eigen values and minimum speed for second fixed
% point (x2)
% this script calculate eigen values from the obtained ch. function
% the matlab defined function eig() use floating pointing approximation
% which provide wrong result sometimes e.g our case in J3

clear all;
clc;
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
g2 = function_g(1/alpha); % g(1/alpha)
%%
c_old = 2*sqrt(mu0*D);  % intial c

temp = 0;
count = 1;
while(temp==0)
    
a1 = c_old + c_old/D;
a2 = ((c_old^2)/D) - g2 - alpha - (mu0/D)*(1 - (1/alpha));
a3 = (c_old/D)*((mu0/alpha) - mu0 - alpha);
a4 = (mu0/D)*(alpha-1);

%the charactersitic equation
p = [1 a1 a2 a3 a4];
% the eigen values
eigJ2_old = roots(p);

% condition for real eigen values
temp_1 = (imag(eigJ2_old(1))<1e-25) &&  (imag(eigJ2_old(2))<1e-25) && (imag(eigJ2_old(3))<1e-25) && (imag(eigJ2_old(4))<1e-25);

%condition for negative real eigen values
% temp_2 = (real(eigJ2_old(1))<0) &&  (real(eigJ2_old(2))<0) && (real(eigJ2_old(3))<0) && (real(eigJ2_old(4))<0);

temp = temp_1; % keep this if only real is enough

% temp = temp_1 && temp_2;  % comment above temp2 and uncomment this for
% real and negative

c_old = c_old + 1e-5;  % increment of c

% j2eigen(count,:) = eigJ2_old;  % uncomment this if eigen values at every
% itreation is required
count = count+1;


end

