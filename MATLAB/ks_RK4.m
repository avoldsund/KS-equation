% Fourth order Runge-Kutta
clear all
close all

f = @(x) cos(x/16).*(1+sin(x/16));

global M k h N
M = 128;
h = (32*pi)/(M);
k = 0.001;
x = 0:h:32*pi;
N = 10000;

U(:,1) = f(x(2:end));


for i = 1:N+1
    
    rk = rkIncrements(U(:,i));
    
    U(:,i+1) = U(:,i) + k/6*(rk(1) + 2*rk(2) + 2*rk(3) + rk(4));
    
end



contourf(U')



figure
mesh(U)