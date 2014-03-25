function feuler_hardcode

f = @(x) cos(x/16) .* (1 + sin(x/16));

L = 32*pi;
M = 2^7;
h = L/M;
x = h:h:L;

h = L/(M-1);

N = 20000;
k = 0.01;

U = zeros(M,N);
U(:,1) = f(x');


    e = ones(M,1);
    diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
    A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

    diagVecB = [-M+1 -1:1 M-1];
    B = (k/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

    diagVecD = [-M+1 -1 1 M-1];
    D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);
    
for n = 1:N
    U(:,n+1) = (eye(M) - A - B)*U(:,n) - D*(U(:,n).^2);
end

figure
contourf(U')


end
