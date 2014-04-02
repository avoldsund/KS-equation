function [] = IMEX_plot()
profile ON
tic
% Function to plot different plots

% The function ref_sol needs the global variables to work
global Ms hs
f = @(x) cos(x/16) .* (1 + sin(x/16));
L = 32*pi;
%f = @(x) (1/sqrt(2)) * sin(x) - (1/8)*sin(2*x);
%L = 2*pi;


M = 2^8;
h = L/M;
x = 0:h:L-h;

N = 2^12;
T = 100;
k = T/N;
t = k*(0:N);

%Reference solution: Must be included for error plots!
Ms = 2^8;  %number of points in reference sol.
hs = L/Ms;
x_ref = 0:hs:L-hs;
U_ref = ref_sol(k,T,x_ref);


%U = spalloc(M,N,M*N);
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
%  titleStr = sprintf('Numerical solution, h = %0.3f, k = %0.3f',h,k)
%  title(titleStr)
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Overview plot of the solution
%  figure
%  [U_c, t_c, x_c] = compress(t, x, U, 1000, 1000);
%  mesh(x_c, t_c, U_c')
%  axis([0 L 0 T])
%  view(-10,45)
%  xlabel('Space [0 32*pi]')
%  ylabelStr = sprintf('Time = %d', T);
%  ylabel(ylabelStr)
%  titleStr = sprintf('Numerical solution, h = %0.3f, k = %0.3f',h,k);
%  title(titleStr)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Error plot between the reference solution and our numerical approximation
%  figure
%  error = zeros(M,N); 
%  error = abs(U_ref(1:Ms/M:Ms,:) - U);
%  [err_c, t_c, x_c] = compress(t, x, error, 1000, 1000);
%  mesh(x_c, t_c, err_c')
%  titleStr = sprintf('Error plot: approximate solution (h = %0.3f) vs reference solution (h = %0.3f). k = %0.3f.',h,hs,k);
%  title(titleStr)
%  ylabelStr = sprintf('Time = %d',T)
%  ylabel(ylabelStr)
%  xlabel('Space [0, 32*pi]')
%  axis([0 L 0 T])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Comparison of the reference solution and the numerical approx. at time T
% figure
% plot(x, U_ref(1:Ms/M:Ms,end),'r', x, U(:,end),'linewidth', 2.5)
% str = sprintf('Comparison: Numerical solution (h = %0.3f) vs reference solution (h = %0.3f) at time %d. k = %0.3f.',h,hs,T,k);
% title(str);
% xlabel('Space [0, 32*pi]')
% ylabelStr = sprintf('U(x,%d)',T);
% ylabel(ylabelStr)
% axis([0 32*pi -4 4])
% legend('Reference solution', 'Numerical solution')
profile OFF
profile viewer
toc
end