% Rk4


close all
clear all


f = @(x) cos(x/16).*(1+sin(x/16));


M = 10;
h = (32*pi)/(M-1);
k = 0.0001;

x = 0:(M-1)*h;
N = 1000;


t = linspace(0,(N+1)*k,N+1);
U(:,1) = f(x);

e = ones(M,1);

s = [-1:1 M-1];
A = spdiags([e, -2*e, e, e], s,M,M);
p = [M-2 M-1 1
g = @(U) 