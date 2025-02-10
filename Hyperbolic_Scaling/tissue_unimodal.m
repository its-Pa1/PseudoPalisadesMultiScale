function q  = tissue_unimodal(x,y,theta,delta, kappa, center_x_DW, center_y_DW)
% This function compute the Non-symmetric fiber distribution(directed mesoscopic tissue (q)), given by eq.
% 2.5.5 of "Diss_Kumar_Pawan.pdf" present in the parent directory

% the d(x,y)
d = 0.25*(exp(-0.005*(x-center_x_DW).*(x-center_x_DW)))-0.25*(exp(-0.005*(y-center_y_DW).*(y-center_y_DW)));

% as the D_w consists only digonal elements, its eigen values are:
lambda_1 = 0.5-d;
lambda_2 = 0.5+d;

% Fraction anisotropy of D_w
FA =  abs(lambda_1-lambda_2)/sqrt((lambda_1^2)+ (lambda_2^2));
% k and N2 as defined in the paper, kappa controls the FA of D_w
k= kappa*FA;
N2 = 1./((2*pi)*besseli(0,k));
% phii computes the leading eigen value
if (d<=0)
    phii = [1,0];
else
    phii = [0,1];
end

% for first term in q for hyperpolic scaling
% first unimodal von mises ditribution
k1 = 0.05*exp(-0.000001*(((x-center_x_DW).^2)+((y-center_y_DW).^2)));
gamma = [1/sqrt(2);1/sqrt(2)];
N22 = 1./((2*pi)*besseli(0,k1));

q2 = N22*exp(k1*[cos(theta),sin(theta)]*gamma);

% finally q, here delta combine both uniomodals von mises
% distribution
q = delta*q2+ (1-delta)*N2.*((exp(k.*phii*[cos(theta);sin(theta)])));


end