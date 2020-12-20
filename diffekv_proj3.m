%% 1.1 
%Time from t to tend
M = 200;
t = 0;
tend = 0.1;
dt = (tend - t)/M;
tt = linspace (0,tend,M+1);

%space from 0 to 1 N points
%larger N => CFL condition fails, CFL 0.1-0.2 ish
N = 30;
xx = linspace(0,1,N+2)';
dx = 1/(N+1);
x = xx(2:end-1);

%initial condition function
%g = @(x) x.^3 + 1;
g = @(x) -x.*(x-1);

%Toeplitz
%Tdx = (diag(-2*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1))/dx^2;
e = ones(N, 1);
T = spdiags([e, -2*e, e], [-1, 0, 1], N, N);
Tdx = T/dx^2;

init = g(x);

unew = zeros(N, M+1);
unew(:,1) = init;


for i = 1:M
    unew(:,i+1) = eulerstep(Tdx, unew(:,i), dt);
end

unew = [zeros(1,M+1); unew; zeros(1,M+1)];
[T,X] = meshgrid(tt,xx);

figure;
surf(T,X,unew)
shading interp
CFL = dt/dx^2

%% 1.2
%Time from t to tend
M = 70;
t = 0;
tend = 0.3;
dt = (tend - t)/M;
tt = linspace (0,tend,M+1);

%space from 0 to 1 N points
N = 300;
xx = linspace(0,1,N+2)';
x = xx(2:end-1);
dx = 1/(N+1);
CFL = dt/dx^2;

%initial condition function
g = @(x) x.^3.*(-x.^2+1);

%Toeplitz
e = ones(N, 1);
T = spdiags([e, -2*e, e], [-1, 0, 1], N, N);

Tdx = T/dx^2;

init = g(x);

unew = zeros(N, M+1);
unew(:,1) = init;


for i = 1:M
    unew(:,i+1) = TRstep(Tdx, unew(:,i), dt);
end

unew = [zeros(1,M+1); unew; zeros(1,M+1)];
[T,X] = meshgrid(tt,xx);

surf(T,X,unew)
CFL = dt/dx^2
shading interp


%% 2
clear
clc
%initial condition function
g = @(x) exp(-100*(x-0.5).^2);

%Time from t to tend, M+1 points
M = 50;
t = 0;
tend = 5;
dt = (tend - t)/M;
tt = linspace (0,tend,M+1);

%space from 0 to (N-1)dx, N points
N = 450;
xx = linspace(0,1,N+1)';
x = xx(1:end-1);
dx = 1/(N);

a = 0.02;

CFL = a*dt/dx


init = g(x);
unew = zeros(N,M);
unew(:,1) = init;

for i = 1:M
    unew(:,i+1) = LaxWen(unew(:,i), a*dt/dx);
end

[T,X] = meshgrid(tt,x);


figure
surf(T,X,unew)

figure
plot(tt, rms(unew))

%% 3.1
clear
clc
%Time from t to tend, M+1 points
M = 30;
t = 0;
tend = 1;
dt = (tend - t)/M;
tt = linspace (0,tend,M+1);

%space from 0 to (N-1)dx, N points
N = 40;
xx = linspace(0,1,N+1)';
x = xx(1:end-1);
dx = 1/(N);

%initial condition function
g = @(x) exp(-25*(x-0.5).^2);
%g = @(x) abs(sin(3.*pi.*x));

d = 0.1;
a = 1;
Pe = abs(a/d)
mesh = Pe * dx

init = g(x);
unew = zeros(N,M);
unew(:,1) = init;

for i = 1:M
    %unew(:,i+1) = convdif(unew(:,i), a, d, dt, dx);
    temp = convdif(unew(:,i), a, d, dt, dx);
    temp(1) = temp(end);
    unew(:, i+1) = temp;
end

unew(end+1,:) = unew(1,:);
[T,X] = meshgrid(tt,xx);

figure
surf(T,X,unew)
shading interp

%% 4.1

clear
clc
%Time from t to tend, M+1 points
M = 1500;
t = 0;
tend = 5;
dt = (tend - t)/(M+1);
tt = linspace (0,tend,M+1);

d = 0.01;

%space from 0 to (N-1)dx, N points
N = 200;
xx = linspace(0,1,N+1)';
x = xx(1:end-1);
dx = 1/(N);

%initial condition function
g = @(x) exp(-25.*(x-0.5).^2);
%g = @(x) abs(sin(3.*pi.*x));

%a = 1;
%Pe = abs(a/d)
%mesh = Pe * dx

init = g(xx);
unew(:,1) = init;

for j = 1:M
    temp = visBur(unew(2:end,j), dt, dx, d);
    unew(:, j+1) = [temp(end); temp];
end

[T,X] = meshgrid(tt,xx);

figure
surf(X,T,unew)
shading interp

function unew = visBur(uold, dt, dx, d)
    N = length(uold);
    I = eye(N);
    
    T = diag(-2*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1);
    T(1, N) = 1;
    T(N, 1) = 1;
    Tdx = T/dx^2;
    
    unew = (I-d*dt/2*Tdx)\(LW(uold, dt, dx) + d*dt/2*Tdx*uold);
end

function unew = LW(uold, dt, dx)
    N = length(uold);
    
    Td = diag(ones(1,N-1),1) - diag(ones(1,N-1),-1);
    Td(1, N) = -1;
    Td(N, 1) = 1;
    Td = Td ./(2*dx);
    
    Tsd = diag(-2*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1);
    Tsd(1, N) = 1;
    Tsd(N, 1) = 1;
    Tsd = Tsd /dx^2;
    
    deri = Td;
    secderi = Tsd;

    unew = uold - dt .* uold .* deri * uold + dt^2/2 .* (2*uold.* (deri*uold).^2  + uold.^2.* (secderi * uold));
end

function unew = convdif(u, a, d, dt, dx)
    N = length(u);
    I = eye(N);
    S = diag(1*ones(1,N-1),1) - diag(1*ones(1,N-1),-1);
    S(1, N) = -1;
    S(N, 1) = 1;
    
    T = diag(-2*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1);
    T(1, N) = 1;
    T(N, 1) = 1;
    
    Tdx = d*T/dx^2 - a*S/2/dx;
    
    unew = TRstep(Tdx, u, dt);
end


%for task 2.1
% u^j_n+1 = a*mu/2(1+a*mu)
function unew = LaxWen(u, amu)
    amuh = amu/2;
    N = length(u);
    Tdxx = diag((1-amu^2)*ones(1,N)) + (amuh*(amu+1))*diag(1*ones(1,N-1),1) + (amuh*(amu-1))*diag(1*ones(1,N-1),-1);
    Tdxx(1, N) = amuh*(amu+1);
    Tdxx(N, 1) = amuh*(amu-1);
    
    unew = Tdxx * u;
end


% for task 1.2
% U_n+1 = (I - dt/2 * T_dx) \ (I + dt/2 * T_dx) * U_n
function unew = TRstep(Tdx, uold, dt)
    tdhalf = dt/2;
    I = eye(length(Tdx));
    divide = (I - tdhalf.*Tdx)\(I + tdhalf.*Tdx);
    unew = divide * uold;
end

% For task 1.1
function unew = eulerstep(Tdx, uold, dt)
    unew = uold + dt.*Tdx*uold;
end