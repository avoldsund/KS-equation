function [] = FE_exp()
% Function to plot different plots

% The function ref_sol needs the global variables to work
global Ms hs
f = @(x) cos(x/16) .* (1 + sin(x/16));
L = 32*pi;

M = 2^8;
h = L/M;
x = 0:h:L-h;

N = 2^13;
T = 100;
k = T/N;
t = k*(0:N);

% Reference solution:
Ms = 2^10;  %number of points in reference sol.
hs = L/Ms;
x_ref = 0:hs:L-hs;
U_ref = ref_sol(k,T,x_ref);


U = zeros(M,N);
U(:,1) = f(x);

% IMEX-scheme
A = (k/(h^4))*second_order_matrix(M)*second_order_matrix(M); % Fourth order matrix
B = (k/(h^2))*second_order_matrix(M);
D = (1/2)* k/(2*h)*first_order_central_matrix(M);

F = (speye(M)+A/2+B/2);
G = (speye(M)-A/2-B/2);
    
for n = 1:N
    U(:,n+1) = F\G*(U(:,n)) - F\D*(U(:,n).^2);
    %error(:,n+1) = abs(U_ref(1:Ms/M:Ms,n+1)-U(:,n+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Contour plot (mesh) of the solution
%  figure
%  [U_c, t_c, x_c] = compress(t, x, U, 1000, 1000);
%  mesh(x_c, t_c, U_c')
%  axis([0 L 0 T])
%  xlabel('Space [0 32*pi]')
%  ylabelStr = sprintf('Time = %d', T)
%  ylabel(ylabelStr)
%  titleStr = sprintf('Numerical solution, M = %d, N = %d',M,N)
%  title(titleStr)
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overview plot of the solution
%  figure
%  [U_c, t_c, x_c] = compress(t, x, U, 500, 500);
%  mesh(x_c, t_c, U_c')
%  axis([0 L 0 T])
%  view(-10,45)
%  xlabel('Space [0 32*pi]')
%  ylabelStr = sprintf('Time = %d', T);
%  ylabel(ylabelStr)
%  titleStr = sprintf('Numerical solution, M = %d, N = %d',M,N);
%  title(titleStr)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Error plot between the reference solution and our numerical approximation
%  figure
%  error = zeros(M,N); 
%  error = abs(U_ref(1:Ms/M:Ms,:) - U);
%  [err_c, t_c, x_c] = compress(t, x, error, 500, 500);
%  mesh(x_c, t_c, err_c')
%  titleStr = sprintf('Error plot: approximate solution (M = %d) vs reference solution (M = %d). N = %d.',M,Ms,N);
%  title(titleStr)
%  ylabelStr = sprintf('Time = %d',T)
%  ylabel(ylabelStr)
%  xlabel('Space [0, 32*pi]')
%  axis([0 L 0 T])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Comparison of the reference solution and the numerical approx. at time T
figure
plot(x, U(:,end), x, U_ref(1:Ms/M:Ms,end), 'r', 'linewidth', 3)
str = sprintf('Comparison: Numerical approx. (M = %d) vs reference solution (M = %d) at time %d. N = %d.',M,Ms,T,N);
title(str);
xlabel('Space [0, 32*pi]')
ylabelStr = sprintf('U(x,%d)',T);
ylabel(ylabelStr)
axis([0 32*pi -4 4])
legend('Numerical approx.','Reference solution')

end