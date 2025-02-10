function y = function_g_new(S, S_max,lambda0, lambda1,kp,km)
% this function computes the g(S(x)) involved in the eq. 2.4.16 of 
% "Diss_Kumar_Pawan.pdf" present in the parent directory 

kd = km/kp;
y = ((lambda1*kd)/S_max)/(((kp*(S/S_max))+km+lambda0)*(((S/S_max)+kd)^2));

end