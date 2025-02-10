%%
clear all;
close all;
clc;

% This script calls all the related funtions & solves the non-dim 1D model 2.5.2 of
% "Diss_Kumar_Pawan.pdf" present in the parent directory.
% The corresponding results are included in the Figure 2.16 of the thesis
% In this case the diff. coeff. is given by eq. 2.5.14 of the thesis

tic
%% data (dimensionalised parameters)
scale = 60*60*24*30; % scale month
alpha_old = 1e-11*scale;
beta_old = 1e-9*scale;
s = (10/3600)*scale;
lambda0_old = (0.1)*scale;
lambda1_old = 0.2*scale;
kp_old = 0.004*scale;
km_old = 0.01*scale;
kd_old = km_old/kp_old;
M_max = 0.8;
S_max = 10^(-6.4);
mu0_old = (0.2/(60*60*24))*scale;
delta = 0.2;
kappa = 3;
center_x_DW = 450;
DT_scale = ((s^2)/lambda0_old);
% Ds_old = 25*DT_scale;
Ds_old = 50*DT_scale;
%% non_dim parameters
beta = 1;
alpha = (alpha_old*S_max)/(beta_old*M_max);
mu0 = (mu0_old*S_max)/(beta_old*M_max);
Ds = 1;
lambda0 = lambda0_old/(kp_old*S_max);
lambda1 = lambda1_old/(kp_old*S_max);
kd = kd_old/S_max;
x_scale = sqrt((beta_old*M_max)/(Ds_old*S_max)); %non-dim factor of X
x_end = 1000*x_scale; % length corresponding to 1000
t_end = 10*(beta_old*M_max)/S_max; % time corresponding to 10 months
%% mesh generation, discretization and memory allocation
x = linspace(0,x_end,201);
h = x(2)-x(1);
Lx = length(x);
dt = 10; % looks crazy but doesn't violate the cfl & also saves memory
t = 0:dt:t_end;
M = zeros(Lx,length(t)); % saves M at each time step
S = zeros(Lx,length(t)); % saves S at each time step
flux = zeros(Lx,1);

%% Setting up both the diff matrices
% tumor diff matrix with 1D DT

[D,A] = set_diff_mat_D1(x,dt,DT_scale/(2*Ds_old), x_scale);

% acidity diff matrix
A2  = diag((1/dt) + (2*Ds/(h*h)) * ones(Lx-2,1), 0) + diag(-Ds/(h*h)* ones(Lx-3,1), -1) ...
    + diag(- (Ds/(h*h))* ones(Lx-3,1), 1);
A2(1,1) = A2(1,1)-Ds/(h*h); % BC
A2(end,end) = A2(end,end)-Ds/(h*h); %BC

B_M = zeros(1,size(A,2)); % RHS for tumor
B_S = zeros(1,size(A,2)); % RHS for acidity
%% Initial conditions
radius = 50*x_scale;
sigma  = 25*x_scale;

centerx = 500*x_scale;

IC=zeros(size(x));
IC2=zeros(size(x));

for i=1:length(x)
    
    dx1 =   (x(i)-centerx).^2;
    distance = sqrt(dx1);
    
    if (distance<radius)
        f_val = (0.005/M_max)*exp(-(dx1)/(2*sigma*sigma));% initial tumor density
    else
        f_val = 0;
    end
    IC(i)=f_val;
end

% for acidity
radius = 30*x_scale;
sigma  = 15*x_scale;

for i=1:length(x)
    
    dx1 =   (x(i)-centerx).^2;
    distance = sqrt(dx1);
    
    if (distance<radius)
        f_val = exp(-(dx1)/(2*sigma*sigma));% initial conc of acidity
    else
        f_val = 0;
    end
    
    IC2(i)=f_val;
end

M_old = IC;
S_old = IC2;
%% Time loop
for j= 1:length(t)
    
    % flux
    for i = 2: Lx-1
        flux(i) = (D(i)*function_g_non2(S_old(i),lambda0, lambda1,kd)* (S_old(i+1)-S_old(i))*(1/h))...
            +((D(i)-D(i-1))/h);
    end
    flux(1) = flux(2);
    flux(Lx) = flux(Lx-1);
    
    for i = 2:Lx-1
        % computation of RHS for tumor eq(source+ advection)
        if (flux(i)<=0)
            B_M(i-1) = (M_old(i)/dt) +(1/h)*(flux(i)*M_old(i)-flux(i-1)*M_old(i-1))...% taxis with upwind
                + ((mu0)*(1-M_old(i))*(1-S_old(i))*M_old(i));%source
        else
            B_M(i-1) = (M_old(i)/dt) + (1/h)*(flux(i+1)*M_old(i+1)-flux(i)*M_old(i))...
                +( (mu0)*(1-M_old(i))*(1-S_old(i))*M_old(i));
        end
        
        % RHS for acidity (source only)
        B_S(i-1) = (S_old(i)/dt)+(beta*M_old(i)-alpha*S_old(i));
        
    end
    
    % solutions at interior nodes
    sol_M = (A)\B_M';
    sol_S = (A2)\B_S';
    
    % boundary conditions
    M(2:end-1,j) = sol_M;
    S(2:end-1,j) = sol_S;
    M(1,j) = M(2,j);
    M(Lx,j) = M(Lx-1,j);
    S(1,j) = S(2,j);
    S(Lx,j) = S(Lx-1,j);
    
    % resetting new solution as old for time marching
    M_old = M(:,j);
    S_old = S(:,j);
    
end
%% Plots
% plot of M and S with time and space

% to check wheather Plots_patterns folder exists otherwise it makes a folder Plots_patterns
if not(isfolder('Plots_patterns'))
    mkdir('Plots_patterns')
end


figure(1)
surf(x,t,M')
view(0,90)
colorbar
caxis([0 2e-5]) % as matlab colormap has only few colors, a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 2e-5)
% another option would be to use log colorbar, however that doesn't work
% properly in all cases
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Glioma cells' , 'Fontsize', 15);
axis([250 500 0, t(end)]);
saveas(gcf,'Plots_patterns/Tumor_D_sin2.png');

figure(2)
surf(x,t,S')
view(0,90)
colorbar
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Acidity' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_patterns/Acid_D_sin2.png');

figure(3)
plot(x,D, 'Linewidth', 1.5)
xlabel('X' , 'Fontsize', 15);
ylabel('D_T' , 'Fontsize', 15);
axis([0 x(end) 0 max(D)]);
title('D_T ', 'Fontsize', 15);
saveas(gcf,'Plots_patterns/Diff_D_sin2.png');

toc
