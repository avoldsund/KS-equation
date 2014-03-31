function U_sol = get_reference(k,T)

global Ms hs

Ms = 2^9;
L = 32*pi;
hs = L/Ms;
x = hs * (0:Ms-1);
% time = 10;
% k = time/2^20;

f = @(x) cos(x/16).*(1+sin(x/16));

y0 = f(x);
options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4);
[~,yy] = ode15s('funcKS', 0:k:T, y0, options);
U_sol = yy';

end