function yy = ref_sol(k,T,x)

f = @(x) cos(x/16).*(1+sin(x/16));

y0 = f(x);
options = odeset('AbsTol', 1e-3, 'RelTol', 1e-3);
[~,yy] = ode15s('funcKS', 0:k:T, y0, options);
yy = yy';

end