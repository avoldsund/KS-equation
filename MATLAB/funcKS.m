function F = funcKS(~,y)
    global M h
%     e = ones(M+1,1);
%     
%     s = [-M -1:1 M];
%     A = spdiags([e, e, -2*e, e, e], s,M+1,M+1);
%     p = [-M -M+1 -2:2 M-1 M];
%     B = spdiags([-4*e, e, e, -4*e, 6*e, -4*e, e, e, -4*e],p,M+1,M+1);
%     q = [-M -1 1 M];
%     C = spdiags([e,-e,e,-e],q,M+1,M+1);
% 
%     F = - (1/h^2*A*y + 1/h^4*B*y + 1/(4*h)*C*y.^2);
    
    


    e = ones(M,1);
    diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
    A = (1/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

    diagVecB = [-M+1 -1:1 M-1];
    B = (1/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

    diagVecD = [-M+1 -1 1 M-1];
    C = (1/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);


    F = - (A*y + B*y + C*y.^2);

end