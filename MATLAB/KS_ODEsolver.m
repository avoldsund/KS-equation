function KS_ODEsolver
M = 19;
h = 1/(M+1);
xinit = 32*pi*(1:M)'/M;

ainit = cos(xinit/16).*(1+sin(xinit/16));
y0 = [ainit, xinit];
tspan = [0 1];

sol = ode15s(@F, tspan, y0);
end

function F=funcHePer(t,y)
%
%
global m
mm=m+1;
e=ones(mm,1);
A=mm^2*spdiags([e -2*e e],-1:1,mm,mm);
A(1,mm)=mm^2;
A(mm,1)=mm^2;
%
F=A*y;
end