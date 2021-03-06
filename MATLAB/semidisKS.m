function F = semidisKS(y)
% Semi-Discretization of Kuramoto-Sivashinsky
%   Gives out the derivative f(y,t)
    global M h
    e = ones(M,1);
    
    s = [-M+1 -1:1 M-1];
    A = spdiags([e, e, -2*e, e, e], s,M,M);
    p = [-M+1 -M+2 -2:2 M-2 M-1];
    B = spdiags([-4*e, e, e, -4*e, 6*e, -4*e, e, e, -4*e],p,M,M);
    q = [-M+1 -1 1 M-1];
    C = spdiags([e,-e,e,-e],q,M,M);

    F = - (1/h^2*A*y + 1/h^4*B*y + 1/(4*h)*C*y.^2);

end

