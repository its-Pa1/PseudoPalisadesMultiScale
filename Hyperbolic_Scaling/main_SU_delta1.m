clear all;
close all;
clc;
tic
%This is the main script, where all related functions have been called
%% Grid generation and memory allocations
h = 5; % spatial spacing
x = 0:h:1000; % x domain
y = x; % y,same as x for a square domain
[X,Y] = meshgrid(x,y); % square grid
dt = 0.0005; % time spacing
t_end = 10; % end time
Lx = length(x); % length of x and y
N = (Lx-2)^2; % number of unknowns in the domain
flux_x = zeros(Lx,Lx); % flux in x direction
flux_y = zeros(Lx,Lx); % flux in y direction
advection_x = zeros(Lx-2,Lx-2); % advection term x-component
advection_y = zeros(Lx-2,Lx-2); % advection term y-component
% these three lines to save 10 values of solutions for 1 time interval(month/week/day)
temp = 1/(10*dt);
M_total = zeros(Lx,Lx,t_end*10+1);
S_total = zeros(Lx,Lx,t_end*10+1);
%% data
delta = 1; % delta for small q computation, the parameter to combine two distributions
kappa = 3; % it controls the FA of D_w %
center_x_DW = 450; % position of x for D_w cross
center_y_DW = 450; % position of y for D_w cross

% parameters in the model (all are in \sec time scale then scaled to month)
scale = 60*60*24*30; % scale month
alpha = 1e-11*scale;
beta  = 1e-9*scale;
s = (10/3600)*scale;
lambda0 = 0.1*scale;
lambda1 = 0.2*scale;
kp = 0.004*scale;
km = 0.01*scale;
DT_scale = ((s^2)/lambda0);
D_s =  DT_scale*50;
M_max = 0.8;
S_max = 10^(-6.4);
mu0 = (0.2/(60*60*24))*scale;
epsilon = 1e-5;

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
M_total (:,:,1) = M_old;
S_total (:,:,1) = S_old;

%% calling both the diffusion matrix functions, Q, Eq, DivDT,DivEq
[q,E_q,div_Eq, Q_value, a,b,c,FA,DivDT_x,DivDT_y,A] = set_tumor_diff(x,y,dt, delta, kappa, center_x_DW,center_y_DW,DT_scale,epsilon);
A_s = set_acidity_diff(x,D_s,dt);
Q = reshape(Q_value(2:end-1,2:end-1),N,1);
Ex = reshape(E_q(1,:,:),size(x,2),size(y,2));
Ey = reshape(E_q(2,:,:),size(x,2),size(y,2));

%% Time loop

count = 2;
count2 = 1;
test_ep = 1e-15;

for j = 0:dt:t_end
    
    for ii = 2:length(x)-1
        for jj = 2:length(y)-1
            
            % flux calculation in both directions
            flux_x(ii,jj) = (s*Ex(ii,jj)...
                + epsilon*DT_scale*Ex(ii,jj)*div_Eq(ii-1,jj-1)...
                - DivDT_x(ii-1,jj-1)...
                - function_g_new(S_old(ii,jj),S_max,lambda0, lambda1,kp,km)*...
                ((a(ii,jj)*((S_old(ii+1,jj)-S_old(ii,jj))/(h))) + ...
                (b(ii,jj)*((S_old(ii,jj+1)-S_old(ii,jj))/(h)))));
            
            flux_y(ii,jj) = (s*Ey(ii,jj)...
                + epsilon*DT_scale*Ey(ii,jj)*div_Eq(ii-1,jj-1) + ...
                - DivDT_y(ii-1,jj-1)...
                - function_g_new(S_old(ii,jj),S_max,lambda0, lambda1,kp,km)*...
                ((b(ii,jj)*((S_old(ii+1,jj)-S_old(ii,jj))/(h))) + ...
                (c(ii,jj)*((S_old(ii,jj+1)-S_old(ii,jj))/(h)))));
        end
    end
    flux_x(1,:)   = flux_x(2,:);          flux_y(1,:)   = flux_y(2,:);
    flux_x(end,:) = flux_x(end-1,:);      flux_y(end,:) = flux_y(end-1,:);
    flux_x(:,1)   = flux_x(:,2);          flux_y(:,1)   = flux_y(:,2);
    flux_x(:,end) = flux_x(:,end-1);      flux_y(:,end) = flux_y(:,end-1);
    
    for kk = 2:Lx-1
        for ll = 2:Lx-1
            % upwinding in both directions (2nd order upwind with van-Leer flux limitation)
            if ((kk==2)||(kk==Lx-1)||(ll==2)||(ll==Lx-1))
                if(flux_x(kk,ll)>=0)
                    advection_x(kk-1,ll-1) = (flux_x(kk,ll)*M_old(kk,ll) - flux_x(kk-1,ll)*M_old(kk-1,ll))/h;
                else
                    advection_x(kk-1,ll-1) = (flux_x(kk+1,ll)*M_old(kk+1,ll) - flux_x(kk,ll)*M_old(kk,ll))/h;
                end
                
                if(flux_y(kk,ll)>0)
                    advection_y(kk-1,ll-1) = (flux_y(kk,ll)*M_old(kk,ll) - flux_y(kk,ll-1)*M_old(kk,ll-1))/h;
                else
                    advection_y(kk-1,ll-1) = (flux_y(kk,ll+1)*M_old(kk,ll+1) - flux_y(kk,ll)*M_old(kk,ll))/h;
                end
            else
                rx_m = (M_old(kk+1,ll)- M_old(kk,ll) + test_ep)/(M_old(kk,ll) - M_old(kk-1,ll) + test_ep);
                rx_p = (M_old(kk,ll) - M_old(kk+1,ll) + test_ep)/(M_old(kk+1,ll)- M_old(kk+2,ll) + test_ep);
                ry_m = (M_old(kk,ll+1)- M_old(kk,ll) + test_ep)/(M_old(kk,ll) - M_old(kk,ll-1) + test_ep);
                ry_p = (M_old(kk,ll) - M_old(kk,ll+1) + test_ep)/(M_old(kk,ll+1)- M_old(kk,ll+2) + test_ep);
                phi_x_m = (abs(rx_m)+rx_m)/(1+abs(rx_m));
                phi_x_p = (abs(rx_p)+rx_p)/(1+abs(rx_p));
                phi_y_m = (abs(ry_m)+ry_m)/(1+abs(ry_m));
                phi_y_p = (abs(ry_p)+ry_p)/(1+abs(ry_p));
                
                if(flux_x(kk,ll)>=0)
                    advection_x(kk-1,ll-1) = ((flux_x(kk,ll)* (M_old(kk,ll)+0.5*phi_x_m*(M_old(kk,ll)-M_old(kk-1,ll))))/h)...
                        - ((flux_x(kk-1,ll)* (M_old(kk-1,ll)+0.5*phi_x_m*(M_old(kk-1,ll)-M_old(kk-2,ll))))/h);
                else
                    advection_x(kk-1,ll-1) = ((flux_x(kk+1,ll)*(M_old(kk+1,ll)-0.5*phi_x_p*(M_old(kk+2,ll)-M_old(kk+1,ll))))/h)...
                        - (flux_x(kk,ll)*(M_old(kk,ll)-0.5*phi_x_p*(M_old(kk+1,ll)-M_old(kk,ll))))/h;
                end
                
                if(flux_y(kk,ll)>=0)
                    advection_y(kk-1,ll-1) = ((flux_y(kk,ll)* (M_old(kk,ll)+0.5*phi_y_m*(M_old(kk,ll)-M_old(kk,ll-1))))/h)...
                        - ((flux_y(kk,ll-1)* (M_old(kk,ll-1)+0.5*phi_y_m*(M_old(kk,ll-1)-M_old(kk,ll-2))))/h);
                else
                    advection_y(kk-1,ll-1) = ((flux_y(kk,ll+1)*(M_old(kk,ll+1)-0.5*phi_y_p*(M_old(kk,ll+2)-M_old(kk,ll+1))))/h)...
                        - (flux_y(kk,ll)*(M_old(kk,ll)-0.5*phi_y_p*(M_old(kk,ll+1)-M_old(kk,ll))))/h;
                end
            end
        end
    end
    B_M = reshape(M_old(2:end-1,2:end-1),N,1);
    B_S = reshape(S_old(2:end-1,2:end-1),N,1);
    advection =  reshape(advection_x,N,1) + reshape(advection_y,N,1);
    
    source_M = epsilon*(mu0)*((1-(B_M/M_max)).*(1-(B_S/S_max))).*B_M; % source glioma cells
    
    
    RHS_M = B_M - dt*advection+ dt*source_M; % entire RHS for tumor eq(source and advection) for IMEX method
    
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
    
    % % %     plots per iteration
    % % %     uncomment till title, if plots at each iteration is required (for debug)
    % % %         figure(1)
    % % %         set(gcf, 'Position',  [100, 600, 500, 500])
    % % %         surf(x,y,M_old')
    % % %         view(0,90)
    % % %         colorbar
    % % %         shading interp
    % % %         colormap jet
    % % %         drawnow
    % % %         title(['Tumor at time step ',num2str(j)], 'Fontsize', 15);
    % % %         figure(2)
    % % %         set(gcf, 'Position',  [600, 600, 500, 500])
    % % %         surf(x,y,S_old')
    % % %         view(0,90)
    % % %         colorbar
    % % %         shading interp
    % % %         colormap jet
    % % %         drawnow
    % % %         title(['Acidity at time step ',num2str(j)], 'Fontsize', 15);
    %
    if (mod(count2,temp)==0)
        M_total (:,:,count) = M_old;
        S_total (:,:,count) = S_old;
        count = count+1;
        
    end
    
    S_old = S_new; % resetting new solution as old for time marching
    M_old = M_new;
    
    
    count2 = count2+1;
end

%% videos
% These two loops save the solution video per each 3 days

% In this section Videos are saved in Videos_delta1 folder
% To check if Videos_delta1 folder exists otherwise to create
if not(isfolder('Videos_delta1'))
    mkdir('Videos_delta1')
end

% to check which video profile supports available in the machine
% if mp4 is not supported then avi format will be used
profiles = VideoWriter.getProfiles();
check_mp4_support = find(ismember({profiles.Name},'MPEG-4'));

if isempty(check_mp4_support)
    video_ext = '.avi';
    v_pro  = 'Motion JPEG AVI';
else
    video_ext = '.mp4';
    v_pro = 'MPEG-4';
end

videofile = VideoWriter(strcat('Videos_delta1/Tumor_hyper', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(M_total,3)
    figure(3)
    %     set(gcf, 'Position',  [100, 600, 500, 500])
    surf(x,y,M_total(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial glioma cells', 'Fontsize', 15);
    else
        title(['Glioma cells at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 15);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F);
end
close(videofile);

videofile = VideoWriter(strcat('Videos_delta1/Acidity_hyper', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(S_total,3)
    figure(4)
    %     set(gcf, 'Position',  [600, 600, 500, 500])
    surf(x,y,S_total(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial acidity', 'Fontsize', 15);
    else
        title(['Acidity at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 15);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);

%% This loop saves the results for each 30 days in folder Plots_delta1 in png
% or eps format with days in file names

% to check wheather Plots folder exists otherwise it makes a folder Plots
if not(isfolder('Plots_delta1'))
    mkdir('Plots_delta1')
end

for i = 1:10:size(M_total,3)
    figure(5)
    %     set(gcf, 'Position',  [100, 600, 500, 500])
    surf(x,y,M_total(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial glioma cells', 'Fontsize', 15);
    else
        title(['Glioma cells at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 15);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots_delta1/Tumor_HS_delta_1_%ddays',3*(i-1)),'epsc');%eps
    saveas(gcf,sprintf('Plots_delta1/Tumor_HS_delta_1_%ddays.png',3*(i-1)));%png
    
    figure(6)
    %     set(gcf, 'Position',  [600, 600, 500, 500])
    surf(x,y,S_total(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial acidity', 'Fontsize', 15);
    else
        title(['Acidity at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 15);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots_delta1/Acidity_HS_delta_1_%ddays',3*(i-1)),'epsc');
    saveas(gcf,sprintf('Plots_delta1/Acidity_HS_delta_1_%ddays.png',3*(i-1)))
end
%% Other plots (FA and Q)

figure(7)
% set(gcf, 'Position',  [100, 600, 500, 500])
surf(x,y,FA','LineStyle','none');
view(0,90);
colorbar
shading interp
colormap('jet')
title('FA')
xlabel('X' , 'Fontsize', 15);
ylabel('Y' , 'Fontsize', 15);
saveas(gcf,'Plots_delta1/FA_HS_delta_1.png');
%saveas(gcf,'Plots_delta1/FA_HS_delta_1','epsc');


figure(8)
% set(gcf, 'Position',  [100, 600, 500, 500])
surf(x,y,Q_value','LineStyle','none');
view(0,90);
colorbar
shading interp
colormap('jet')
title('Q')
xlabel('X' , 'Fontsize', 15);
ylabel('Y' , 'Fontsize', 15);
saveas(gcf,'Plots_delta1/Q_macro_HS_delta_1.png');
% saveas(gcf,'Plots_delta1/Q_macro_HS_delta_1','epsc');

figure(9)
% set(gcf, 'Position',  [100, 600, 500, 500])
surf(x,y,q(:,:,26)','LineStyle','none');
view(0,90);
colorbar
shading interp
colormap('jet')
title('$q_h$ at $\xi = \frac{\pi}{2}$','Fontsize', 18,'Interpreter','latex');
xlabel('X' , 'Fontsize', 15);
ylabel('Y' , 'Fontsize', 15);
saveas(gcf,'Plots_delta1/q_HS_delta_1.png');
% saveas(gcf,'Plots_delta1/q_HS_delta_1','epsc');

figure(10)
quiver(X,Y,Ex',Ey',1.5);
axis([0,1000,0,1000])
title('E_q')
xlabel('X' , 'Fontsize', 15);
ylabel('Y' , 'Fontsize', 15);
saveas(gcf,'Plots_delta1/Eq_HS_delta_1.png');

figure(11)
mask1 = (x <=550) & (x>=400);
mask2 = (y <=550) & (y>=400);
quiver(X(mask1,mask2),Y(mask1, mask2),Ex(mask1,mask2)',Ey(mask1,mask2)',2);
axis([400,550,400,550])
title('Zoomed E_q')
xlabel('X' , 'Fontsize', 15);
ylabel('Y' , 'Fontsize', 15);
saveas(gcf,'Plots_delta1/Eq_zoomed_HS_delta_1.png');

%% uncomment to save the workspace
if not(isfolder('Mat_files'))
    mkdir('Mat_files')
end
save('Mat_files/main_2D_HS_delta_1.mat');
%%
toc