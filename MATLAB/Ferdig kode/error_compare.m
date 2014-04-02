function [e_e e_i N_j] = error_compare


f = @(x) cos(x/16) .* (1 + sin(x/16));
L = 32*pi;

 global Ms hs
%   %number of points in reference sol.
% 


M = 2^7;
h = L/M;
k = 0.01;
Ms = 2^10;
hs = L/Ms;
y = 0:hs:L-hs;
iter = 20;

N_j = zeros(iter,1);
e_e = zeros(iter,1);
e_i = zeros(iter,1);

for j = 1:iter
    N = 1000*j;

    x = 0:h:L-h;
    T = k*N;
    
    yy = ref_sol(k,T,y);

    U_i = zeros(M,N);
    U_i(:,1) = f(x);
    U_e = U_i;

    %Generating matrices

    A = k/(h^2)*second_order_matrix(M);
    B = k/(h^4)*second_order_matrix(M)*second_order_matrix(M);
    D = k/(4*h)*first_order_central_matrix(M);

    F = (speye(M)+0.5*A+0.5*B);
    G = (speye(M)-0.5*A-0.5*B);


    % Time step N iterations:
    for n = 1:N-1
         U_e(:,n+1) = (eye(M)-A-B)*U_e(:,n) - 0.5*D*(U_e(:,n).^2);
    end
    e_e(j) = norm(yy(1:Ms/M:Ms,N)-U_e(:,N), Inf);


    % Time step N iterations:
    for n = 1:N-1
        U_i(:,n+1) = F\G*(U_i(:,n)) - F\D*(U_i(:,n).^2);
    end
    e_i(j) = norm(yy(1:Ms/M:Ms,N)-U_i(:,N), Inf);
    
    N_j(j) = N;
end

figure
plot(N_j, e_i, 'ro-',  N_j, e_e,'bo-');
legend('Implicit scheme', 'Explicit scheme')


end

