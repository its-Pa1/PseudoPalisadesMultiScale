function [DC,A] = set_diff_mat_alex(x,dt, DT_scale, x_scale)
%%
% This function deals with assembling the diffusion matrix required for IMEX discretization of
% 1D tumor eq of the model 2.5.2 of
% "Diss_Kumar_Pawan.pdf" present in the parent directory
% In this case the diffusion coefficeint is D_2(x) given by 2.5.15


% outputs:
% DC = 1D DT
% A = tumor diffusion matrix

% inputs:
% x: the space variable
% dt = time spacing
% DT_scale = s^2/lambda0
% x_scale = non-dim scale
h = x(2)-x(1);
DC = zeros(size(x,2),1);
Lx = length(x);
dia_alpha = zeros(Lx-2,1);
dia_beta = zeros(Lx-2,1);
a_d = 1; % aplha for the d(x)

%% Computation of DT
for i = 1:length(x)
    if ((250*x_scale<=x(i)) && (x(i)<=500*x_scale))
        DC(i) = (sin(0.04*pi*(x(i)-250*x_scale)))^(2*a_d);
    elseif ((750*x_scale<=x(i)) && (x(i)<=1000*x_scale))
        DC(i) = (sin(0.04*pi*(x(i)-750*x_scale)))^(2*a_d);
    end
end

DC = DT_scale*DC;
%% tumor diff matrix
for i = 2:Lx-1
    dia_alpha(i-1) = (1/(h*h))*0.5*(DC(i)+ DC(i+1));
    dia_beta(i-1) = (1/(h*h))*0.5*(DC(i)+ DC(i-1));
end
diag_1 = dia_beta;
diag_1(1:end-1) = dia_beta(2:end);
diag0 = (dia_alpha + dia_beta);
diag0(1) = diag0(1)-dia_beta(1);
diag0(end) = diag0(end)-dia_alpha(end);
diag1 = dia_alpha;
diag1(2:end) = diag1(1:end-1);
A = spdiags([-diag_1,(1/dt)+diag0, -diag1],[-1,0,1],Lx-2,Lx-2);
% A1 = full(A);
end