tic
% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time
% U(n+1) = (I - A - B)U(n) - D*R(U(n))

f = @(x) cos(x/16) .* (1 + sin(x/16));

% Initial values
L = 32*pi;
M = 2^7;
h = L/(M);
x = 0:h:32*pi-h;

N = 5000;
k = 0.01;
T = 5;
% N_time = T/k;

r = k/(h^4)

% Creating the U-matrix and inserting boundary condition

U = zeros(M, N);
U(:,1) = f(x);
U_R = U;

% Construction of the A-, B- and D-matrix
e = ones(M,1);
diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

diagVecB = [-M+1 -1:1 M-1];
B = (k/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

diagVecD_R = [-M+1 -1 1 M-1];
D = (k/(2*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD_R, M, M);

% D = spdiags([1*e -1*e 1*e -1*e], diagVecD_R, M, M)*D;
% diagVecD = [-M+1 -1 1 M-1];
% D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);


% spectral_radius = max(abs(eig(C1)));



% Testing out: U(x,0) = -2*U_xx - U_xxxx
% rho = 1/(2^7)*cos(x/16).*(1 + 4*sin(x/16)) - 1/(2^16)*cos(x/16).*(1+16*sin(x/16));
% R = spdiags(rho', 0, M, M);
% 


% Finding U(T=20)
% for n = 1:N-1
%     U(:,n+1) = (eye(M) - A - B)*U(:,n) - D*(U(:,n).^2);
% end

% R_diag = (U(:,1) + U(:,N_time))/2;
R_diag = (U(:,1));
R = spdiags(R_diag, 0, M, M);
% 
C = (eye(M) - A - B); 

E = full(D*R);
% plot(eig(E), '*');

% C1 = (eye(M)-A-B-D);
% C2 = (eye(M)-A-B);
eigen_C = (eig(C));
% plot(eig(C2), '*')
% spec_C2 = max(abs(eig(C2)));
% 
% 
for n = 1:N-1
    U_R(:,n+1) = (eye(M) - A - B)*U_R(:,n) - 0.5*D*R*U_R(:,n);
end

% 
% C2 = (eye(M) - A - B - D*R);
% spec_C2 = max(abs(eig(C2)));
% err = (U-U_R);
% norm(err)

% figure(1)
% contourf(U')
% figure(2)
contourf(U_R')

toc