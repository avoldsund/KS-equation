function FE_exp

global Ms hs
f = @(x) cos(x/16) .* (1 + sin(x/16));
L = 32*pi;

Ms = 2^8;  %number of points in reference sol.
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

B = second_order_matrix(M,h);
A = k*B*B;          % Fourth order matrix
B = k*B;
D = k*first_order_central_matrix(M,h);
    
for n = 1:N-1
    U(:,n+1) = (eye(M) - A - B)*U(:,n) - D*(U(:,n).^2);
    error(:,n+1) = abs(yy(1:Ms/M:Ms,n+1)-U(:,n+1));
end
% figure(1)
% err = abs(yy(:,N)-U(:,N));
% plot(err, 'r');
% hold on
% plot(yy(:,N), 'b');
% hold on
 %plot(U(:,N), 'g');
norm(error)
figure(2)
%meshc(error')

% figure
 contourf(U')


end