function err = FE_imp

global Ms hs
f = @(x) cos(x/16) .* (1 + sin(x/16));
L = 32*pi;

Ms = 2^10;  %number of points in reference sol.
hs = L/Ms;
y = 0:hs:L-hs;

M = 2^7;
h = L/M;
x = 0:h:L-h;

N = 20000;
k = 0.01;
T = k*N;

yy = ref_sol(k,T,y);

U = zeros(M,N);
error = zeros(M,N);
U(:,1) = f(x);

%Generating matrices for imp.exp. method

A = k/(2*h^2)*second_order_matrix(M);

B = k/(2*h^4)*second_order_matrix(M)*second_order_matrix(M);

D = k/(4*h)*first_order_central_matrix(M);

F = (speye(M)+A+B);
G = (speye(M)-A-B);


% Time step N iterations:
for n = 1:N-1
    U(:,n+1) = F\G*(U(:,n)) - F\D*(U(:,n).^2);
    %error(:,n+1) = abs(yy(1:Ms/M:Ms,n+1)-U(:,n+1));
end

%err = norm(error);
err = norm(yy(1:Ms/M:Ms,N)-U(:,N));
%contourf(U')


end