function [q,Q_value, a,b,c,FA,DivDT_x,DivDT_y,A] = ...
    set_tumor_diff(x,y,dt, delta,kappa, center_x_DW,center_y_DW, DT_scale)
% This function deals with everything required for tumor diffusion
% computation
% the non-negative discretization (the stencils are shwon in Table 2.3  of 
% "Diss_Kumar_Pawan.pdf" present in the parent directory) has been used for 
% tumor diff tensor

% outputs:
% q = mesoscopic tissue(q)
% Q_value = macroscopic tissue(Q)
% a = DT(1,1) first element of first row of tumor diff tensor
% b = DT(1,2) = DT(2,1) 
% c = DT(2,2)
% FA = fractional anisotropy
% DivDT_x = x component of div of DT
% DivDT_y = y component of div of DT
% A = tumor diffusion matrix

% inputs:
% x,y : the space variables
% dt = time spacing
% delta = required for computing small q, it combines the two distributions
% kappa = it controls the FA of D_w
% center_x_DW = x position of the cross for D_w(as our D_w makes a cross type sign)
% center_y_DW = y position of the cross for D_w
% DT_scale = s^2/lambda0

%% memory allocations
Lx = length(x);
N = (Lx-2)^2;
h = x(2)-x(1);
% theta = linspace(0,pi,100);
theta = linspace(0,pi,101);
q = zeros(size(x,2),size(y,2),size(theta,2));
Q_value = zeros(size(x,2),size(y,2));
DC = zeros(2,2,size(x,2),size(y,2));
FA = zeros(size(x,2), size(y,2));
alpha1 = zeros(N,1);
alpha2 = zeros(N,1);
alpha3 = zeros(N,1);
alpha4 = zeros(N,1);
alpha5 = zeros(N,1);
alpha6 = zeros(N,1);
alpha7 = zeros(N,1);
alpha8 = zeros(N,1);
alpha9 = zeros(N,1);
a = zeros(size(x,2),size(y,2));
b = zeros(size(x,2),size(y,2));
c = zeros(size(x,2),size(y,2));
DivDT_x = zeros(Lx-2,Lx-2);
DivDT_y = zeros(Lx-2,Lx-2);

%% Loops to save the Q and q at all points and theta(small q depends of theta)
for j = 1:length(y)
    for i = 1:length(x)
        Q_value(i,j) = tissue_Q_macro(x(i),y(j),center_x_DW,center_y_DW);
        for l = 1:length(theta)
            q(i,j,l) = tissue_q_un(x(i),y(j),theta(l),delta, kappa, center_x_DW, center_y_DW);
        end
        
    end
end
%% Computation of DT and storing its components in a, b and c
% here trapezoidal method is used for numerical integration
for j = 1:length(y)
    for i = 1:length(x)
        sum2 = 0;
        
        for l = 2:length(theta)-1
            sum2 = sum2 + q(i,j,l)*([cos(theta(l));sin(theta(l))]*[cos(theta(l));sin(theta(l))]');
        end
        DC(:,:,i,j) = DT_scale*2*0.5*(theta(2)-theta(1))*(2*sum2 ...
            + q(i,j,1)*[cos(theta(1));sin(theta(1))]*[cos(theta(1));sin(theta(1))]'...
            + q(i,j,end)*[cos(theta(end));sin(theta(end))]*[cos(theta(end));sin(theta(end))]');
        
        temp_d = DC(:,:,i,j);
        a(i,j) = temp_d(1,1);
        b(i,j) = temp_d(1,2);
        c(i,j) = temp_d(2,2);
    end
end
%% Here FA is stored at all points in the domain
for j = 1:length(y)
    for i = 1:length(x)
        lambda = eig(DC(:,:,i,j));
        FA(i,j) = (abs(lambda(1)-lambda(2)))/(sqrt((lambda(1)^2)+(lambda(2)^2)));
    end
end
%% In this section, all the 9 stencils(as mentioned in Weikart et al, pdf in this folder)
%have been saved as alpha's. alpha1 is the upper-left stencil, alpha2:
%upper-middle and so on. Also, in these loops, the div of DT has been
%stored

for j = 2:Lx-1
    for i = 2:Lx-1
        k = (j-2)*(Lx-2)+i-1;
        ii = i-1;
        jj = j-1;
        
        alpha1(k) = ((abs(b(i-1,j+1)) - b(i-1,j+1)) + (abs(b(i,j)) - b(i,j)))/(4*h*h);
        
        alpha2(k) = ((c(i,j+1)+c(i,j))/(2*h*h)) - ((abs(b(i,j+1))+abs(b(i,j)))/(2*h*h));
        
        alpha3(k) = ((abs(b(i+1,j+1)) +  b(i+1,j+1))/(4*h*h)) + (abs(b(i,j))+b(i,j))/(4*h*h);
        
        alpha4(k) = ((a(i-1,j)+a(i,j))/(2*h*h)) - ((abs(b(i-1,j))+abs(b(i,j)))/(2*h*h));
        
        alpha5(k) =  -((a(i-1,j)+2*a(i,j)+a(i+1,j))/(2*h*h)) - ((abs(b(i-1,j+1))-b(i-1,j+1)+...
            abs(b(i+1,j+1))+b(i+1,j+1))/(4*h*h)) - ((abs(b(i-1,j-1))+b(i-1,j-1)+...
            abs(b(i+1,j-1))-b(i+1,j-1))/(4*h*h)) + ((abs(b(i-1,j))+abs(b(i+1,j))+...
            abs(b(i,j-1))+abs(b(i,j+1))+2*abs(b(i,j)))/(2*h*h)) - ((c(i,j-1)+2*c(i,j)+...
            c(i,j+1))/(2*h*h));
        
        alpha6(k) = ((a(i+1,j)+a(i,j))/(2*h*h)) - ((abs(b(i+1,j))+abs(b(i,j)))/(2*h*h));
        
        alpha7(k) = ((abs(b(i-1,j-1))+b(i-1,j-1))/(4*h*h)) + (abs(b(i,j))+b(i,j))/(4*h*h);
        
        alpha8(k) = ((c(i,j-1)+c(i,j))/(2*h*h)) - ((abs(b(i,j-1))+abs(b(i,j)))/(2*h*h));
        
        alpha9(k) = ((abs(b(i+1,j-1))-b(i+1,j-1))/(4*h*h)) + ((abs(b(i,j)) - b(i,j))/(4*h*h));
        %         alpha1(k) = i;
        
        DivDT_x(ii,jj) = (a(i,j) - a(i-1,j) + b(i,j) - b(i-1,j))/h;
        DivDT_y(ii,jj) = (b(i,j) - b(i,j-1) + c(i,j) - c(i,j-1))/h;
    end
end



%% As the 9 stencils results in a sparse matrix with 9 digonals
% the diagonal numbering has been done from left to right
% these look crazy, but needed to be done to store in the sparse
% matrix format and to use the matlab function spdiags

% upper diagonals
diag9 = alpha3;
diag9(Lx:end) = diag9(1:end-(Lx-1));
diag9(1:Lx-2:end)=0;

diag8 = alpha2;
diag8(1:Lx-2:end) = alpha2(1:Lx-2:end)+alpha1(1:Lx-2:end);
diag8(Lx-2:Lx-2:end) = alpha2(Lx-2:Lx-2:end)+alpha3(Lx-2:Lx-2:end);
diag8(Lx-1:end) = diag8(1:end-(Lx-2));

diag7 = alpha1;
diag7(Lx-2:end) = diag7(1:end-(Lx-3));
diag7(Lx-2:Lx-2:end) = 0;

diag6 = alpha6;
diag6(1:Lx-2) = alpha6(1:Lx-2)+alpha9(1:Lx-2);
diag6(end-(Lx-2):end) = alpha6(end-(Lx-2):end)+alpha3(end-(Lx-2):end);
diag6(Lx-2:Lx-2:end)=0;
diag6(2:end) = diag6(1:end-1);


%main diagonal
diag5 = alpha5;
diag5(1) = alpha4(1) + alpha5(1)+ alpha7(1)+alpha8(1);
diag5(2:Lx-3) = alpha5(2:Lx-3)+alpha8(2:Lx-3);
diag5(Lx-2) = alpha5(Lx-2)+alpha6(Lx-2)+alpha8(Lx-2)+alpha9(Lx-2);
diag5(Lx-1:Lx-2:N-2*(Lx-2)+1) = alpha4(Lx-1:Lx-2:N-2*(Lx-2)+1)+alpha5(Lx-1:Lx-2:N-2*(Lx-2)+1);
diag5(2*(Lx-2):Lx-2:(Lx-2-1)*(Lx-2)) = alpha5(2*(Lx-2):Lx-2:(Lx-2-1)*(Lx-2))+...
    alpha6(2*(Lx-2):Lx-2:(Lx-2-1)*(Lx-2));
diag5(N-(Lx-2)+1) =  alpha1(N-(Lx-2)+1)+alpha2(N-(Lx-2)+1)+alpha4(N-(Lx-2)+1)+...
    alpha5(N-(Lx-2)+1);
diag5(N-(Lx-2)+2:N-1) = alpha2((N-(Lx-2)+2:N-1)) + alpha5((N-(Lx-2)+2:N-1));
diag5(N) = alpha2(N) + alpha3(N) + alpha5(N) + alpha6(N);


%Lower diagonals

diag1 = alpha7;
diag1(1:end-(Lx-1)) = diag1(Lx:end);
diag1(Lx-2:Lx-2:end) = 0;

diag2 = alpha8;
diag2(1:Lx-2:end) = alpha7(1:Lx-2:end)+alpha8(1:Lx-2:end);
diag2(Lx-2:Lx-2:end) = alpha8(Lx-2:Lx-2:end)+alpha9(Lx-2:Lx-2:end);
diag2(1:end-(Lx-2)) = diag2(Lx-1:end);

diag3 = alpha9;
diag3(1:end-(Lx-3)) = diag3(Lx-2:end);
diag3(1:Lx-2:end) = 0;

diag4 = alpha4;
diag4(1:Lx-2) = alpha4(1:Lx-2)+alpha7(1:Lx-2);
diag4(end-(Lx-3):end) = alpha1(end-(Lx-3):end)+alpha4(end-(Lx-3):end);
diag4(1:end-1) = diag4(2:end);
diag4(Lx-2:Lx-2:end) = 0;

%% finally the diffusion matrix
A = spdiags([-dt*diag1, -dt*diag2, -dt*diag3 , -dt*diag4, (1-dt*diag5), -dt*diag6, ...
    -dt*diag7, -dt*diag8, -dt*diag9],[-(Lx-2+1),-(Lx-2),-(Lx-2-1),-1,0,1,Lx-2-1, Lx-2, Lx-2+1],N,N);

end