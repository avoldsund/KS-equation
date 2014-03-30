%  MOVIE
clear all
close all
global Ms hs
f = @(x) cos(x/16) .* (1 + sin(x/16));
L = 32*pi;

Ms = 2^8;  %number of points in reference sol.
hs = L/Ms;
y = 0:hs:L-hs;

N = 10000;
k = 0.01;
T = k*N;


y0 = f(y);
options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[~,yy] = ode15s('funcKS', 0:k:T, y0, options);

figure
for i=1:size(yy,1)
    plot(y,[yy(i,end);yy(i,:)'])
    drawnow
end

figure
mesh(y,tt,[yy(:,end) yy])
xlabel('rom')
ylabel('tid')