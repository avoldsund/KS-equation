function [U U_ref error] = kuramotoSivashinsky(method, errPlot)
%% KURAMOTO-SIVASHINSKY(method, errorPlot)
%   A function for solving the Kuramoto-Sivashinsky equation with periodic
%   boundary conditions: U(0,t) = U(L,t), and initial conditions: U(x,0) = f(x).
%   The initial conditions used are f(x) = cos(x/16) * (1 + sin(x/16)) on the
%   interval L = 32*pi.
%   INPUT: Method of choice. 1 for explicit method (FE-type) and 2 for
%   implicit-explicit method (CN-type). ErrorPlot: 1 for a plot comparing
%   the numerical solution with the reference solution.
%   OUTPUT: The numerical solution (U), the reference solution (U_ref) and
%   the error matrix (error)
%   NOTE: To change the length of time steps (k), or the length of the
%   steps in space (h) for either the reference solution or the numerical
%   solution, it is possible to change them in KS_solver. NB: To maintain
%   stability for the explicit method, a requirement for k/h^4 < 1/8 is
%   necessary.

% Solving the KS-equation
[U U_ref x t M M_ref] = KS_solver(method);

% Plotting the solution. The function uses a compress function to save time
plotKS(t, x, U, 500, 500, method)

% Error plot
if errPlot == 1
    % Because both M_ref and M are in the form 2^k, we are able to compare
    % the matrices, even though they are different sizes.
    error = abs((U_ref(1:M_ref/M:M_ref, :)) - U);
    
    % Plots an error plot
    errorPlot(error, t, x)
end
end

function [U U_ref x t M M_ref] = KS_solver(method)
%% KS_SOLVER(method)
%   Solves the KS-equation using an implicit-explicit scheme.
%   OUTPUT: [U U_ref x x_ref t error] where U is the numerical solution. U_ref is the
%   reference solution aquired by using ode15s.

% Initial conditions and interval length
f = @(x) cos(x/16) .* (1 + sin(x/16));
L = 32*pi;

% M: Number of space steps
% h: Step length in space
% x: The x-grid
M = 2^7;
h = L/M;
x = 0:h:L-h;

% N: Number of time steps
% T: The time where the KS-equation is solved
% k: Length of time step
% t: The time grid
T = 100;
k = 0.01;
N = T/k;
t = k * (0:N-1);

% Variables for the reference solution
M_ref = 2^10;  %number of points in reference sol.
h_ref = L/M_ref;
x_ref = 0:h_ref:L-h_ref;

% Calculating the reference solution
U_ref = ref_sol(k, T, x_ref, M_ref, h_ref);

U = zeros(M,N);
U(:,1) = f(x);

% IMEX-scheme
A = (k/(h^4))*second_order_matrix(M)*second_order_matrix(M); % Fourth order matrix
B = (k/(h^2))*second_order_matrix(M);
D = (1/2)* k/(2*h)*first_order_central_matrix(M);

% Iteration for the time steps, different solutions depending on method
if method == 1
    for n = 1:N-1
        U(:,n+1) = (eye(M)-A-B)*U(:,n) - D*(U(:,n).^2);
    end
else
    F = (speye(M)+A/2+B/2);
    G = (speye(M)-A/2-B/2);
    for n = 1:N-1
        U(:,n+1) = F\G*(U(:,n)) - F\D*(U(:,n).^2);
    end
end
end

function [] = plotKS(t, x, U, xRes, yRes, method)
%% PLOTKS(t, x, U, xRes, yRes, method)
%   plotKS gives a mesh of the solution U. To save time, it compresses the
%   matrix, and instead plots a compressed solution to save time

% The compress function
[U_comp t_comp x_comp] = compress(t, x, U, 500, 500);

% Mesh plot of U
mesh(x_comp, t_comp, U_comp')
xlabel('Space [0, 32*pi]')
ylabel('Time [0, 100]')
axis([0 32*pi 0 100])

if method == 1
    title('Numerical solution of the KS-equation by explicit method')
else
    title('Numerical solution of the KS-equation by implicit-explicit method')
end

end

function [] = errorPlot(error, t, x)
    %Error plot between the reference solution and our numerical approximation
 figure
 [err_c, t_c, x_c] = compress(t, x, error, 750, 750);
 mesh(x_c, t_c, err_c')
 set(gca, 'clim', [0, 4.3525])
 colorbar
 title('Error plot: Difference between numerical solution and reference solution')
 xlabel('space [0, 32*pi]')
 ylabel('time [0, 100]')
 view(0, 90)
end

function A2 = second_order_matrix(M)
%% SECOND_ORDER_MATRIX(M)
%   One of the matrices
    e = ones(M,1);
    diagVecB = [-M+1 -1:1 M-1];
    A2 = spdiags([e e -2*e e e], diagVecB, M, M);

end

function D = first_order_central_matrix(M)
%% FIRST_ORDER_CENTRAL_MATRIX(M
% One of the matrices
     e = ones(M,1);
     diagVecD = [-M+1 -1 1 M-1];
     D = spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);

end

function yy = ref_sol(k,T,x, M_ref, h_ref)
%% REF_SOL(k, T, x, M_ref, h_ref)
    % Function that solves the KS-equation by semi-discretization. Uses the
    % ode15s-solver.

f = @(x) cos(x/16).*(1+sin(x/16));

y0 = f(x);
options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4);
[~,yy] = ode15s(@(t,y) funcKS(t, y, M_ref, h_ref), [0:k:T-k], y0, options);
yy = yy';
end

function F = funcKS(t, y, M_ref, h_ref)
%% FUNCKS(t, y, M_ref, h_ref)
%   Function used by the ref_sol()-function to compute a reference
%   solution.
%   INPUT: y is the initial condition: y = f(x), t the time steps
    
    e = ones(M_ref,1);
    diagVecA = [-M_ref+1 -M_ref+2 -2:2 M_ref-2 M_ref-1];
    A = (1/(h_ref^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M_ref, M_ref);

    diagVecB = [-M_ref+1 -1:1 M_ref-1];
    B = (1/(h_ref^2)) * spdiags([e e -2*e e e], diagVecB, M_ref, M_ref);

    diagVecD = [-M_ref+1 -1 1 M_ref-1];
    C = (1/(4*h_ref)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M_ref, M_ref);

    F = - (A*y + B*y + C*y.^2);
end

function [ V, x, y ] = compress( X, Y, U, xnum, ynum )
%% COMPRESS(U, xnum ynum) 
%       Compresses the numerical grid matrix U into a matrix V of size 
%       xnum x ynum elements. Assumes X and Y are linearly spaced.
%       [Z, xcomp, tcomp] = compress(x, t, U, nx, nt)

%       x, t, U er numerisk resultat, nx er antall elementer du vil ha på x-aksen langs den komprimerte matrisen, 
%       og nt er antall elementer langs t-aksen i den komprimerte matrisen

x = linspace(min(X), max(X), xnum);
y = linspace(min(Y), max(Y), ynum);

[xq, yq] = meshgrid(x, y);
V = interp2(X, Y, U, xq, yq);

assert(size(V, 1) == ynum);
assert(size(V, 2) == xnum);

% Function provided by Andreas Borgen Longva
end