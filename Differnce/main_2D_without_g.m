
%This is the main script, where all related functions have been called
%% Grid generation and memory allocations
h = 5; % spatial spacing
x = 0:h:1000; % x domain
y = x; % y,same as x for a square domain
[X,Y] = meshgrid(x,y); % square grid
dt = .01; % time spacing
t_end = 10; % end time % was 5
Lx = length(x); % length of x and y
N = (Lx-2)^2; % number of unknowns in the domain
flux_x = zeros(Lx,Lx); % flux in x direction
flux_y = zeros(Lx,Lx); % flux in y direction
advection_x = zeros(Lx-2,Lx-2); % advection term x-component
advection_y = zeros(Lx-2,Lx-2); % advection term y-component
% these three lines to save 10 values of solutions for 1 time interval(month/week/day)
temp = 1/(10*dt);
M_total_0g = zeros(Lx,Lx,t_end*10+1);
S_total_0g = zeros(Lx,Lx,t_end*10+1);
%% data
delta = 0.2; % delta for small q computation, the parameter to combine two distributions
kappa = 3; % it controls the FA of D_w %
center_x_DW = 450; % position of x for D_w cross
center_y_DW = 450; % position of y for D_w cross

% parameters in the model (all are in \sec time scale then scaled to month)
scale = 60*60*24*30; % scale month
alpha = 1e-11*scale;
beta  = 1e-9*scale;
s = (10/3600)*scale;
lambda0 = 0.1*scale;
lambda1 = 0.1*scale;
kp = 0.004*scale;
km = 0.01*scale;
DT_scale = ((s^2)/lambda0);
D_s =  DT_scale*50;
M_max = 0.8;
S_max = 10^(-6.4);
mu0 = (0.2/(60*60*24))*scale;

%% Initial conditions
% for tumor
sigma1  = 25;   sigma2  = 20;
centerx = 500; centerx2 = 600;
centery = 500; centery2 = 500;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
M_old       =  0.005*(exp(-exponent) + exp(-exponent2));
M_old = M_old';

% for acidity
sigma1  = 15;   sigma2  = 10;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
S_old = (10^(-7)*exp(-exponent) + 10^(-7)*exp(-exponent2));
S_old = S_old';

% save the initial distributions
M_total_0g (:,:,1) = M_old;
S_total_0g (:,:,1) = S_old;

%% calling both the diffusion matrix functions
[q,Q_value, a,b,c,FA,DivDT_x,DivDT_y,A] = set_tumor_diff(x,y,dt, delta, kappa, center_x_DW,center_y_DW,DT_scale);
A_s = set_acidity_diff(x,D_s,dt);
Q = reshape(Q_value(2:end-1,2:end-1),N,1);
%% Time loop
count = 2;
count2 = 1;
t = 0:dt:t_end;
flux_x_test = zeros(1,length(t));
flux_y_test = zeros(1,length(t));
for j = 1:length(t)
    
    for ii = 2:length(x)-1
        for jj = 2:length(y)-1
            % flux calculation in both directions
            flux_x(ii,jj) =  DivDT_x(ii-1,jj-1);
            
            flux_y(ii,jj) =  DivDT_y(ii-1,jj-1);
        end
    end
    flux_x(1,:)   = flux_x(2,:);          flux_y(1,:)   = flux_y(2,:);
    flux_x(end,:) = flux_x(end-1,:);      flux_y(end,:) = flux_y(end-1,:);
    flux_x(:,1)   = flux_x(:,2);          flux_y(:,1)   = flux_y(:,2);
    flux_x(:,end) = flux_x(:,end-1);      flux_y(:,end) = flux_y(:,end-1);
    
    flux_x_test(j) = max(max(flux_x));
    flux_y_test(j) = max(max(flux_y));
    
    for kk = 2:Lx-1
        for ll = 2:Lx-1
            % upwinding in both directions
            if(flux_x(kk,ll)<=0)
                advection_x(kk-1,ll-1) = (flux_x(kk,ll)*M_old(kk,ll) - flux_x(kk-1,ll)*M_old(kk-1,ll))/h;
            else
                advection_x(kk-1,ll-1) = (flux_x(kk+1,ll)*M_old(kk+1,ll) - flux_x(kk,ll)*M_old(kk,ll))/h;
            end
            
            if (flux_y(kk,ll)<=0)
                advection_y(kk-1,ll-1) = (flux_y(kk,ll)*M_old(kk,ll) - flux_y(kk,ll-1)*M_old(kk,ll-1))/h;
            else
                advection_y(kk-1,ll-1) = (flux_y(kk,ll+1)*M_old(kk,ll+1) - flux_y(kk,ll)*M_old(kk,ll))/h;
            end
        end
    end
    
    B_M = reshape(M_old(2:end-1,2:end-1),N,1);
    B_S = reshape(S_old(2:end-1,2:end-1),N,1);
    advection =  reshape(advection_x,N,1) + reshape(advection_y,N,1);
    
    source_M = (mu0)*((1-(B_M/M_max)).*(1-(B_S/S_max))).*B_M; % source for glioma cells
    
    
    RHS_M = B_M + dt*advection+ dt*source_M; % entire RHS for tumor eq(source and advection) for IMEX method
    
    RHS_S = B_S+ dt*beta*B_M - dt*alpha*B_S;% RHS containing source & uptake for acidity eq
    
    sol = A\RHS_M; % solution for tumor at interior nodes
    sol2 = A_s\RHS_S; % solution for acidity at interior nodes
    
    % Boundary conditions
    M_sol = reshape(sol,Lx-2,Lx-2);
    new = [M_sol(1,1:end);M_sol;M_sol(end,1:end)];
    M_new = [new(1:end,1),new, new(1:end,end)];
    
    S_sol = reshape(sol2,Lx-2,Lx-2);
    new_S = [S_sol(1,1:end);S_sol;S_sol(end,1:end)];
    S_new = [new_S(1:end,1),new_S, new_S(1:end,end)];
    

    if (mod(count2,temp)==0)
        M_total_0g (:,:,count) = M_old;
        S_total_0g (:,:,count) = S_old;
        count = count+1;
        
    end
    
    S_old = S_new; % resetting new solution as old for time marching
    M_old = M_new;
    
    count2 = count2+1;
end