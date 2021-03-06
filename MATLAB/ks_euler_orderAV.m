function ks_euler_orderAV
% KURAMOTO-SIVASHINSKY-EQUATION WITH FORWARD DIFFERENCE
% ORDER OF CONVERGENCE

% L = 32*pi
f = @(x) cos(x/16).*(1+sin(x/16));
% Reference solution:
it = 6;
M_sol = 2^10+1;
k = 0.0001;
T = 10;
N = floor(T/k);
u_sol = reference_sol(M_sol, N, k);
err_norm = zeros(it,1);
h_p = zeros(it,1);


for j = 1:it
    
    M = 2^j+1;
    x = linspace(0,32*pi,M);
    h = (32*pi)/(M-1);

    
    % Initializing matrices for equation
    % U(n+1) = U(n) - AxxU(n) - AxxxxU(n) - delta(U(n)^2) :
    
    w = ones(1,M-1);
    Axx = sparse(-2*eye(M) + diag(w,1) + diag(w,-1));
    Axx(1,M) = 1;
    Axx(M,1) = 1;
    Axx = k/(h^2)*Axx;
    
    u = ones(1,M-2);
    v = -4*ones(1,M-1);
    Axxxx = sparse(6*eye(M)+diag(v,1)+diag(v,-1)+diag(u,2)+diag(u,-2)+diag([1;1],M-2)+diag([1;1],-M+2));
    Axxxx(1,M) = -4;
    Axxxx(M,1) = -4;
    Axxxx = k/(h^4)*Axxxx;
    
    D = sparse(zeros(M) + diag(w,1) + diag(-w,-1));
    D(M,1) = 1;
    D(1,M) = -1;
    D = k/(4*h)*D;
    
    % Initial value of U:
    U = f(x);
    % Time step N iterations:
    for n = 1:N
        U(n+1,:) = (eye(M) - Axx - Axxxx)*(U(n,:)') - D*(U(n,:)'.^2);%
    end
    mesh(U)
    pause

    error = zeros(1,M);

        for i = 0:M-1
            error(i+1) = U(N,i+1) - u_sol(N,(M_sol-1)/(M-1)*i+1);
        end
    
        
        
        
    err_norm(j) = norm(error, Inf);
    h_p(j) = h;
end

loglog(h_p, err_norm,'r', h_p, h_p)
end


function U = reference_sol(M, N, k)
	f = @(x) cos(x/16).*(1+sin(x/16));
    x = linspace(0,32*pi,M);
    h = (32*pi)/(M-1);
    
    % Initializing matrices for equation
    % U(n+1) = U(n) - AxxU(n) - AxxxxU(n) - delta(U(n)^2) :
    
    w = ones(1,M-1);
    Axx = sparse(-2*eye(M) + diag(w,1) + diag(w,-1));
    Axx(1,M) = 1;
    Axx(M,1) = 1;
    Axx = k/(h^2)*Axx;
    
    u = ones(1,M-2);
    v = -4*ones(1,M-1);
    Axxxx = sparse(6*eye(M)+diag(v,1)+diag(v,-1)+diag(u,2)+diag(u,-2)+diag([1;1],M-2)+diag([1;1],-M+2));
    Axxxx(1,M) = -4;
    Axxxx(M,1) = -4;
    Axxxx = k/(h^4)*Axxxx;
    
    D = sparse(zeros(M) + diag(w,1) + diag(-w,-1));
    D(M,1) = 1;
    D(1,M) = -1;
    D = k/(4*h)*D;
    
    % Initial value of U:
    U = f(x);

    % Time step N iterations:
    for n = 1:N
        U(n+1,:) = (eye(M) - Axx - Axxxx)*(U(n,:)') - D*(U(n,:)'.^2);%
    end
    %mesh(U)
    %pause
end
