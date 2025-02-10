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
M_total = zeros(Lx,Lx,t_end*10+1);
S_total = zeros(Lx,Lx,t_end*10+1);
%% data
delta = 0.2; % delta for small q computation, the parameter to combine two distributions
kappa = 3; % it controls the FA of D_w %
center_x_DW = 450; % position of x for D_w cross
center_y_DW = 450; % position of y for D_w cross

% parameters in the model (all are in \sec time scale then scaled to month)
scale = 60*60*24*30; % scale month
alpha = 1e-6*scale;
beta  = 1e-10*scale;
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

%% Initial conditions
% for tumor
sigma1  = 25;   sigma2  = 20;
centerx = 500; centerx2 = 600;
centery = 500; centery2 = 500;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
M_old       =  0.005*(exp(-exponent)+exp(-exponent2));
M_old = M_old';

% for acidity
sigma1  = 15;   sigma2  = 10;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
S_old = (10^(-7)*exp(-exponent)+10^(-7)*exp(-exponent2));
S_old = S_old';

% save the initial distributions
M_total (:,:,1) = M_old;
S_total (:,:,1) = S_old;

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
            flux_x(ii,jj) =  DivDT_x(ii-1,jj-1)+function_g_new(S_old(ii,jj), S_max,lambda0, lambda1,kp,km)*...
                (a(ii,jj)*((S_old(ii,jj)-S_old(ii-1,jj))/h) + ...
                b(ii,jj)*((S_old(ii,jj)-S_old(ii,jj-1))/h));
            
            flux_y(ii,jj) =  DivDT_y(ii-1,jj-1)+function_g_new(S_old(ii,jj), S_max,lambda0, lambda1,kp,km)*...
                (b(ii,jj)*((S_old(ii,jj)-S_old(ii-1,jj))/h) + ...
                c(ii,jj)*((S_old(ii,jj)-S_old(ii,jj-1))/h));
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
    
    source_M = (mu0)*((1-(B_M/M_max)).*(1-(B_S/S_max))).*B_M; % source glioma cells
    
    
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
    
    % plots per iteration
    % uncomment till title, if plots at each iteration is required (for debug)
    % % %    figure(1)
    % % %     set(gcf, 'Position',  [100, 600, 500, 500])
    % % %     surf(x,y,M_old')
    % % %     view(0,90)
    % % %     colorbar
    % % %     shading interp
    % % %     colormap jet
    % % %     drawnow
    % % %     title(['Tumor at time step ',num2str(j)], 'Fontsize', 15);
    % % %     figure(2)
    % % %     set(gcf, 'Position',  [600, 600, 500, 500])
    % % %     surf(x,y,S_old')
    % % %     view(0,90)
    % % %     colorbar
    % % %     shading interp
    % % %     colormap jet
    % % %     drawnow
    % % %     title(['Acidity at time step ',num2str(j)], 'Fontsize', 15);
    
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

% In this section Videos are saved in Videos_lower_grades folder
% To check if Videos_lower_grades folder exists otherwise to create
if not(isfolder('Videos_lower_grades'))
    mkdir('Videos_lower_grades')
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

videofile = VideoWriter(strcat('Videos_lower_grades/Tumor', video_ext),v_pro);
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

videofile = VideoWriter(strcat('Videos_lower_grades/Acidity', video_ext),v_pro);
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

%% This loop saves the results for each 30 days in folder Plots_lower_grades in png
% or eps format with days in file names

% to check wheather Plots folder exists otherwise it makes a folder Plots
if not(isfolder('Plots_lower_grades'))
    mkdir('Plots_lower_grades')
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
    %     saveas(gcf,sprintf('Plots_lower_grades/Tumor_lower_grades_%ddays',3*(i-1)),'epsc');%eps
    saveas(gcf,sprintf('Plots_lower_grades/Tumor_lower_grades_%ddays.png',3*(i-1)));%png
    
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
    %     saveas(gcf,sprintf('Plots_lower_grades/Acidity_lower_grades_%ddays',3*(i-1)),'epsc');
    saveas(gcf,sprintf('Plots_lower_grades/Acidity_lower_grades_%ddays.png',3*(i-1)))
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
saveas(gcf,'Plots_lower_grades/FA_lower_grades.png');
%saveas(gcf,'Plots_lower_grades/FA_lower_grades','epsc');

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
saveas(gcf,'Plots_lower_grades/Q_macro_lower_grades.png');
% saveas(gcf,'Plots_lower_grades/Q_macro_lower_grades','epsc');

figure(9)
% set(gcf, 'Position',  [100, 600, 500, 500])
surf(x,y,q(:,:,51)','LineStyle','none');
view(0,90);
colorbar
shading interp
colormap('jet')
title('q at $\xi = \frac{\pi}{2}$','Fontsize', 18,'Interpreter','latex');
xlabel('X' , 'Fontsize', 15);
ylabel('Y' , 'Fontsize', 15);
saveas(gcf,'Plots_lower_grades/q_lower_grades.png');
% saveas(gcf,'Plots_lower_grades/q_lower_grades','epsc');

%% uncomment to save the workspace
if not(isfolder('Mat_files'))
    mkdir('Mat_files')
end
save('Mat_files/main_2D_lower_grades.mat');
%%
toc