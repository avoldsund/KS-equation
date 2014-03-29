function D = first_order_central_matrix(M)

     e = ones(M,1);
     diagVecD = [-M+1 -1 1 M-1];
     D = spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);

end