function Q_value  = tissue_Q_macro(x,y, center_x_DW, center_y_DW)
% this function computes the macroscopic tissue density(Q) as described on
% page 33 of section 2.5 of 
% "Diss_Kumar_Pawan.pdf" present in the parent directory 

% d(x,y) for D_w
d = 0.25*(exp(-0.005*(x-center_x_DW).*(x-center_x_DW)))-0.25*(exp(-0.005*(y-center_y_DW).*(y-center_y_DW)));

% leading eigen value
if (d<=0)
    lambda = 0.5-d;
else
    lambda = 0.5+d;
end

% as our D_w is simple, we have an explicit formula for Q
% Q_value = 1 - (1./(8*(lambda.^(3/2)))); in 3D
Q_value = 1 - (1./(4*lambda));% in 2D
end