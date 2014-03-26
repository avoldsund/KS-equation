function F = funcKS(~,y)
    global Ms h
    
    e = ones(Ms,1);
    diagVecA = [-Ms+1 -Ms+2 -2:2 Ms-2 Ms-1];
    A = (1/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, Ms, Ms);

    diagVecB = [-Ms+1 -1:1 Ms-1];
    B = (1/(h^2)) * spdiags([e e -2*e e e], diagVecB, Ms, Ms);

    diagVecD = [-Ms+1 -1 1 Ms-1];
    C = (1/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, Ms, Ms);


    F = - (A*y + B*y + C*y.^2);
 
end