function [error] = FE_exp

global Ms hs
%f = @(x) cos(x/16) .* (1 + sin(x/16));
f = @(x) (1/sqrt(2))*sin(x) -(1/8)*sin(2*x)
L = 32*pi
%L = 32*pi;

Ms = 2^10;  %number of points in reference sol.
hs = L/Ms;
y = 0:hs:L-hs;

M = 2^6;
h = L/M
x = 0:h:L-h;

N = 2^18;
T = 3.05;
k = T/N
t = k*(0:N-1);
r = (k/h^4)

% yy = ref_sol(k,T,y);

U = zeros(M,N);
% error = zeros(M,N);
U(:,1) = f(x);

B = second_order_matrix(M);
A = (k/h^4)*B*B;          % Fourth order matrix
B = (k/(h^2))*B;
D = k/(2*h)*first_order_central_matrix(M);
    
for n = 1:N-1
    U(:,n+1) = (eye(M) - A - B)*U(:,n) - 0.5*D*(U(:,n).^2);
%     error(:,n+1) = abs(yy(1:Ms/M:Ms,n+1)-U(:,n+1));
end

figure
[UNew, xNew, yNew] = compress(t, x, U, 500, 500);
 mesh(UNew)
disp('lol')
% figure(1)
%  err = norm(yy(1:Ms/M:Ms,N)-U(:,N));
% plot(err, 'r');
% hold on
% plot(yy(:,N), 'b');
% hold on
 %plot(U(:,N), 'g');
 
 
 % Finding error plot (difference between ref. and our solution)
 % Using the compress-function
% [errNew, xNew, yNew] = compress(t, x, error, 2000, 2000);

% figure
% mesh(xNew,yNew,errNew')
% view(0,90) %Contour plot
% %view(90,270) Contour plot of the solution
% title('Error plot: approximate solution (M = 256) vs reference solution (M = 1024). N = 2^{18}.')
% ylabel('time')
% xlabel('space [0, 32*pi]')
% axis([0 32*pi 0 100])
 
% err = norm(error);
% size(error)
% mesh(x,t,error')
% xlabel('space')
% ylabel('time')
% view(0,90)
% contourf(error')

end