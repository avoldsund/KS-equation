
function error_norm = space_convergence
% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time
% Reference solution with ode15


f = @(x) cos(x/16).*(1+sin(x/16));

global Ms hs
L = 32*pi;
Ms = 2^10;
hs = L/(Ms);
k = 0.0001;
x0 = 0:hs:L-hs;
N = 500000;
T = N*k;
%size(x0)
y0 = f(x0);


options=odeset('AbsTol',1e-3,'RelTol',1e-3);
[tt,yy] = ode15s('funcKS',0:k:T,y0,options);
%[space, time] = size(yy');

%figure
%contourf(yy);
yy = yy';

min = 6;
max = 9;
num = max-min+1;
error_norm = zeros(num,1);
h_p = zeros(num,1);

% figure
% contourf(yy')
% pause
% close all

for j = min:max

    M = 2^j;
    h = L/(M);
    %k/(h^2) <= 0.0026 for convergence
    % x = h * (1:M)
    % x = 32*pi/(M-1) * (1:M)
    x = 0:h:L-h;
    
    
    % Boundary conditions
    % u(x,0) = f(x)
    % u(0,t) = u(L,t)

    f = @(x) cos(x/16) .* (1 + sin(x/16));

    % Creating the U-matrix and inserting boundary condition
    U = zeros(M, N);
    U(:,1) = f(x');

    % Construction of matrices for U(n+1) = (I - A - B)*U(n) - D*(U(n)).^2 
    e = ones(M,1);
    diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
    A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

    diagVecB = [-M+1 -1:1 M-1];
    B = (k/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

    diagVecD = [-M+1 -1 1 M-1];
    D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);

    % Iteration over time
    for n = 1:N
        U(:,n+1) = (eye(M)-A-B)*U(:,n) - D*(U(:,n).^2);
    end

%     figure
%     contourf(U')
%     pause
%     close all

    % error compared to reference solution
    error = zeros(1,M);
    for i = 0:M-1
        error(i+1) = U(i+1,N) - yy((Ms)/(M)*i+1,N);
    end
    error_norm(j-min+1) = norm(error, Inf);
    h_p(j-min+1) = h;
end

figure
    loglog(h_p, error_norm, 'r', h_p, h_p, 'b', h_p, h_p.^2, 'g');
    
end
