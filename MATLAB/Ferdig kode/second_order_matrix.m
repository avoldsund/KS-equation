function A2 = second_order_matrix(M, h)

    e = ones(M,1);
    diagVecB = [-M+1 -1:1 M-1];
    A2 = (1/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

end