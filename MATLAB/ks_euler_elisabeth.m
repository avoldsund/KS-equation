tic
% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time
% U(n+1) = (I - A - B)U(n) - D*(U(n)^2)

% Initial values
L = 32*pi;
M = 128;
h = L/(M);
x = 0:h:L-h;
N = 20000;
k = 0.01;

f = @(x) cos(x/16) .* (1 + sin(x/16));

% Creating the U-matrix and inserting boundary condition
U = zeros(M, N);
U(:,1) = f(x');

% Construction of the A-, B- and D-matrix
    e = ones(M,1);
    diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
    A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

    diagVecB = [-M+1 -1:1 M-1];
    B = (k/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

    diagVecD = [-M+1 -1 1 M-1];
    D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);


% Running 
for n = 1:N
    U(:,n+1) = (eye(M) - A - B)*U(:,n) - D*(U(:,n).^2);
end

figure
%contourf(t, x, U)
contourf(U')
%mesh(U')
toc


%%

% function reference
% 
% f = @(x) cos(x/16).*(1+sin(x/16));
% 
% M = 128;
% h = (32*pi)/(M+1);
% k = 0.01;
% x = 0:h:32*pi;
% N = 20000;
% T = N*k;
% 
% y0 = f(x(2:end));
% 
% options=odeset('AbsTol',1e-3,'RelTol',1e-3);
% [tt,yy] = ode15s('funcKS',[0 T],y0,options);
% 
% figure
% for i=1:size(yy,1)
%     plot(x,[yy(i,end);yy(i,:)'])
%     drawnow
% end
% 
% figure
% mesh(x,tt,[yy(:,end) yy])
% xlabel('rom')
% ylabel('tid')
% 
% figure
% contourf(yy)
% end