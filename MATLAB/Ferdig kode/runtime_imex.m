function [t_e t_i N_j] = runtime_imex


f = @(x) cos(x/16) .* (1 + sin(x/16));
L = 32*pi;

% global Ms hs
% Ms = 2^10;  %number of points in reference sol.
% hs = L/Ms;
% y = 0:hs:L-hs;

M = 2^7;
k = 0.01;

iter = 10;
%T = k*N;
N_j = zeros(iter,1);
t_e = zeros(iter,1);
t_i = zeros(iter,1);

for j = 1:iter
    N = 2000*j;
    h = L/M;
    x = 0:h:L-h;
    % yy = ref_sol(k,T,y);

    U_i = zeros(M,N);
    U_i(:,1) = f(x);
    U_e = U_i;

    %Generating matrices

    A = k/(h^2)*second_order_matrix(M);
    B = k/(h^4)*second_order_matrix(M)*second_order_matrix(M);
    D = k/(4*h)*first_order_central_matrix(M);

    F = (speye(M)+0.5*A+0.5*B);
    G = (speye(M)-0.5*A-0.5*B);


    tic
    % Time step N iterations:
    for n = 1:N-1
         U_e(:,n+1) = (eye(M)-A-B)*U_e(:,n) - 0.5*D*(U_e(:,n).^2);
    end
    t_e(j) = toc;


    tic
    % Time step N iterations:
    for n = 1:N-1
        U_i(:,n+1) = F\G*(U_i(:,n)) - F\D*(U_i(:,n).^2);
    end
    t_i(j) = toc;
    
    N_j(j) = N;
end

figure
plot(N_j, t_i, 'ro-',  N_j, t_e,'bo-');
legend('Implicit scheme', 'Explicit scheme')
% err_i = norm(yy(1:Ms/M:Ms,N)-U_i(:,N), Inf);
% err_e = norm(yy(1:Ms/M:Ms,N)-U_e(:,N), Inf);



end

