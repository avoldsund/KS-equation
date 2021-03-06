function feuler_hardcode

global M h L N k
f = @(x) cos(x/16) .* (1 + sin(x/16));

L = 32*pi;
M = 2^7;
h = L/(M);
x = 0:h:L-h;

N = 50000;
k = 0.01;
T = N*k;

% y0 = f(x);
% options=odeset('AbsTol',1e-3,'RelTol',1e-3);
% [tt,yy] = ode15s('funcKS',[0:k:T],y0,options);
% yy = yy';


U = zeros(M,N);
error = zeros(M,N);
U(:,1) = f(x);

    e = ones(M,1);
    diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
    A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

    diagVecB = [-M+1 -1:1 M-1];
    B = (k/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

    diagVecD = [-M+1 -1 1 M-1];
    D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);
    
    eig(eye(M) - A - B)
    
% for n = 1:N-1
%     U(:,n+1) = (eye(M) - A - B)*U(:,n) - D*(U(:,n).^2);
%     error(:,n+1) = abs(yy(:,n+1)-U(:,n+1));
% end
% figure(1)
% err = abs(yy(:,N)-U(:,N));
% plot(err, 'r');
% hold on
% plot(yy(:,N), 'b');
% hold on
 %plot(U(:,N), 'g');

figure(2)
meshc(error')

% figure
% contourf(U')


end
