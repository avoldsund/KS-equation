% Spatial grid
N = 64; L = 2*pi; x = ([0:N-1]*2*L/N)'; dt = 2^(-10);

% Spectral differentiation matrices
D1 = i*(pi/L)*[0:N/2-1 0 -N/2+1:-1];
D2 = D1.^2; D2((N/2)+1) = -(N*pi/(2*L))^2;
D4 = D2.^2; c = -D2-D4;

% Evaluating the coefficients of the ETD4RK method
% Using Cauchy integral formula
R = 1; N1 = 32; 
r = R*exp(2*i*pi*(1:N1)/N1);
c1 = c*dt; c2 = c1/2; 
E1 = exp(c1); E = exp(c2);

for k = 1:N
C1(k) = real(mean((dt/2)*((exp(c2(k)+r)-1)./(c2(k)+r))));
C2(k) = real(mean(dt*((-4-c1(k)-r+exp(c1(k)+r).*(4-3*(c1(k)+r)+(c1(k)+r).^2))...
./(c1(k)+r).^3)));
C3(k) = real(mean(dt*((2+c1(k)+r+(c1(k)+r-2).*exp(c1(k)+r))./(c1(k)+r).^3)));

C4(k) = real(mean(dt*((-4-3*(c1(k)+r)-(c1(k)+r).^2+(4-c1(k)-r).*exp(c1(k)+r)) ...
./(c1(k)+r).^3)));
end

% Initial condition
u = exp(cos(x/2)); uhat = fft(u);

% Solve PDE:
tmax = 60; nmax = round(tmax/dt); nc = 60; nplt = floor(nmax/nc);
udata = u; tdata = 0;
min1 = min(u); max1 = max(u);
for n = 1:nmax
    t = n*dt;
    uhat1_x = D1.*fft(u.^2)/2;
    ahat = (E.*uhat)-(C1'.*uhat1_x); a=real(ifft(ahat));
    bhat = (E.*uhat)-(C1'.*D1.*fft(a.^2)/2); b=real(ifft(bhat));
    chat = (E.*ahat)-(C1'.*(D1.*fft(b.^2)-uhat1_x)); C=real(ifft(chat));
    uhat = (E1.*uhat)-(C2'.*uhat1_x+C3'.*D1.*(fft(a.^2)+fft(b.^2)) ...
    +C4'.*D1.*fft(C.^2)/2);
    u = real(ifft(uhat));
    if mod(n,nplt) == 0
    udata = [udata u]; 
    tdata = [tdata t];
    min1 = [min1 min(u)]; 
    max1 = [max1 max(u)];
    end
end

% plot results:
set(gcf,'renderer','zbuffer'), clf, drawnow
mesh(x,tdata,udata'), colormap(1e-6*[1 1 1])
xlabel x, ylabel t, zlabel u, grid on
axis([0 2*L 0 tmax floor(min(min1)) ceil(max(max1))])
set(gca,'ztick',[floor(min(min1)) ceil(max(max1))])