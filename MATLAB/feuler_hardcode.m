function feuler_hardcode

global M h L N k
f = @(x) cos(x/16) .* (1 + sin(x/16));

L = 32*pi;
M = 2^7;
h = L/M;
x = h:h:L;

h = L/(M-1);

N = 20000;
k = 0.01;

u = ref_sol;

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

contourf(U')


end


%%
function yy = ref_sol

f = @(x) cos(x/16).*(1+sin(x/16));

global M h L N k

h = L/M;
k = 0.01;
x = h:h:L;
N = 20000;
T = N*k;

y0 = f(x');

options=odeset('AbsTol',1e-3,'RelTol',1e-3);
[~,yy] = ode15s('funcKS',0:k:T,y0,options);

end

function F = funcKS(~,y)
    global M h 
    
    e = ones(M,1);
    diagA = [-M+1 -M+2 -2:2 M-2 M-1];
    A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagA, M, M);

    diagB = [-M+1 -1:1 M-1];
    B = (k/(h^2)) * spdiags([e e -2*e e e], diagB, M, M);

    diagD = [-M+1 -1 1 M-1];
    D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagD, M, M);

    F = - (1/h^2*A*y + 1/h^4*B*y + 1/(4*h)*D*y.^2);
end