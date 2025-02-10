function y = function_g_non2(S,lambda0, lambda1,kd)
% this function computes the non-dimensionalized g(S(x)) as mentioned in subsection 2.5.1 of 
% "Diss_Kumar_Pawan.pdf" present in the parent directory 

y = ((lambda1*kd))./((S+kd+lambda0).*((S+kd).^2));

end