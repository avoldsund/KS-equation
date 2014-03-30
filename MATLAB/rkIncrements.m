function increments = rkIncrements(y)
% Computing Runge Kutta Increments using the function semidisKS
%   k_1 is the increment based on the slope at the beginning of the interval
%   k_2 is the increment based on the slope at the midpoint of the interval
%   k_3 is again the increment based on the slope at the midpoint
%   k_4 is the increment based on the slope at the end of the interval
global k

k_1 = semidisKS(y);

k_2 = semidisKS(y + k/2*k_1);

k_3 = semidisKS(y + k/2*k_2);

k_4 = semidisKS(y + k*k_3);

increments = [k_1 k_2 k_3 k_4];

end

