% Rk4

close all
clear all


f = @(x) cos(x/16).*(1+sin(x/16));

global M k h N
M = 100;
h = (32*pi)/(M);
k = 0.001;
x = 0:h:32*pi;
N = 10000;
T = N*k;

y0 = f(x(2:end));

options=odeset('AbsTol',1e-10,'RelTol',1e-9);
[tt,yy] = ode15s('funcKS',[0 T],y0,options);

%figure
%for i=1:size(yy,1)
 %   plot(x,[yy(i,end);yy(i,:)'])
 %   drawnow
%end

figure
mesh(x,tt,[yy(:,end) yy])
xlabel('rom')
ylabel('tid')

figure
contourf(yy)