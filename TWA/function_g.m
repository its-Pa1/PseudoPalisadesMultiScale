function y = function_g(S)
% this function computes the non-dimensionalized g(S(x)) as mentioned in subsection 2.5.1 of 
% "Diss_Kumar_Pawan.pdf" present in the parent directory 

kp = 0.004;
km = 0.01;
kd = km/kp;
lambda0_old = 0.1;
lambda1_old = 0.01;
lambda0 = lambda0_old/kp;
lambda1 = lambda1_old/kp;

y = lambda1*kd/( ((S + kd)^2)*(S + kd + lambda0));
end