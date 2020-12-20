f = @(t, y) -2*y;
told = 0;
uold = 1;

T = 2; %steget man tar från told
Ns = 2.^(1:10); %alla olika N, antal steg mellan told och T
hs = T./Ns;
exact = exp(-2*T);
error = [];
for k = 1:length(Ns)
    N = Ns(k);
    h = hs(k);
    ys = uold;
    for n = 1:N
        ys(:, n+1) = RK4step(f, told+n*h, ys(:,n), h); 
    end
    error(k) = norm(ys(:,N+1) - exact);
end


loglog(hs, error)
hold on;
grid on;
loglog(hs, hs.^4*error(end)/hs(end).^4)
xlabel('h');
ylabel('Global error');
%% del 2
t0 = 0;
tf = 2;
y0 = 1;
tol = 0.00001;
%[y, error] = RK34step(f, told, uold, h);
[t,y] = adaptiveRK34(f, t0, tf, y0, tol);
figure;
hold on;
grid on;
plot(t,y)
xlabel('t');
ylabel('y');
y2 = exp(-2.*t);
plot(t,y2)


%% del 3
lotka = @(t, u) [3*u(1)-9*u(1)*u(2); 15*u(1)*u(2)-15*u(2)];
H = @(x, y) 15.*x + 9.*y - 15.*log(x) - 3.*log(y);
t0 = 0;
tf = 150;
y0 = [[1; 1]];
[t,y] = adaptiveRK34(lotka, t0, tf, y0, 10^(-11));
figure;
plot(t,y);
xlabel('Time')
ylabel('H(x,y)')

HXY = H(y(1,:), y(2,:));

drift = abs(HXY./H(y(1,1),y(2,1))-1);

figure
loglog(t,drift);
xlabel('Time')
ylabel('Drifting error')


%% del 4.1
vanderpol = @(t, u) [u(2); 100.*(1-u(1).^2).*u(2)-u(1)];
t0 = 0;
tf = 2000;
y0 = [[ 2 ; 0]];
[t,y] = adaptiveRK34(vanderpol, t0, tf, y0, 10^(-8));
figure;
plot(t,y(2,:));
figure;
plot(y(1,:),y(2,:));

%% del 4.2
mys = [10 15 22 33 47 68 100 150 220 330 470 680 1000];


y0 = [[ 2 ; 0]];
N = [];
for i = 1:length(mys)
    my = mys(i);
    vanderpol = @(t, u) [u(2); my.*(1-u(1).^2).*u(2)-u(1)];
    t0 = 0;
    tf = 0.7*my;
    [t,y] = adaptiveRK34(vanderpol, t0, tf, y0, 10^(-8));
    N(i) = length(y);
end

%q ~ 2, stiffness ~ C*my^2
loglog(mys,N);
hold on;
loglog(mys, mys.^2*N(end)/mys(end).^2)

%% 4.3
mys = [10 15 22 33 47 68 100 150 220 330 470 680 1000000];
y0 = [[ 2 ; 0]];
N = [];
for i = 1:length(mys)
    my = mys(i);
    vanderpol = @(t, u) [u(2); my.*(1-u(1).^2).*u(2)-u(1)];
    t0 = 0;
    tf = 0.7*my;
    [t,y] = ode15s(vanderpol, [t0 tf] , y0);
    N(i) = length(y);
end

%q ~ 2, stiffness ~ C*my^2
loglog(mys,N);
hold on;

%it is very quick




function unew = RK4step(f, told, uold, h)
b = [1/6; 1/3; 1/3; 1/6];
A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
c = [0; 0.5; 0.5; 1];

Yprim = [];
% Y'_1 = f(t_n, y_n)
Yprim(:,1) = f(told, uold);
for i = 1:size(A)-1 %Y'_i = f(t_n + h*c_i*Y'_i-1)
    Yprim(:,i+1) = f(told + h*c(i+1), uold + h*c(i+1)*Yprim(:,i));
end

unew = uold + h*Yprim*b;

end

function [unew, err] = RK34step(f, told, uold, h)
b = [1/6; 1/3; 1/3; 1/6];
A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
c = [0; 0.5; 0.5; 1];
cz = [0; 1/2; 1];
bz = [1/6; 2/3; 1/6];
Az = [0 0 0; 0.5 0 0; -1 2 0];

Yprim = [];
% Y'_1 = f(t_n, y_n)
Yprim(:,1) = f(told, uold);
for i = 1:size(A)-1 %Y'_i = f(t_n + h*c_i*Y'_i-1)
    Yprim(:,i+1) = f(told + h*c(i+1), uold + h*c(i+1)*Yprim(:,i));
end
Zprim = f(told + h, uold - h*Yprim(:,1)+2*h*Yprim(:,2));

unew = uold + h*Yprim*b;
err = norm(h/6*(2*Yprim(:,2)+Zprim - 2*Yprim(:,3) - Yprim(:,4)));

end

function hnew = newstep(tol, err, errold, hold, k)
    hnew = ((tol)/(err))^(2/(3*k))*((tol)/(errold))^(-1/(3*k))*hold;
end




function [t,y] = adaptiveRK34(f, t0, tf, y0, tol)

y = [[y0]];
t = [t0];
h = [(abs(tf-t0)*tol^(1/4))/(100*(1+norm(f(t0, y0))))];
i = 1;
err = [[tol]];
while t(end)<tf
    i = i+1;
    [y(:,i), err(:,i)] = RK34step(f, t(i-1), y(:,i-1), h(i-1));
    err(:,i) = norm(err(:,i));
    h(i) = newstep(tol, err(:,i), err(:,i-1), h(i-1), 4);
    t(i) = t(i-1) + h(i-1);
end

h(i) = tf-t(i-1);
t(i) = tf;
[y(:,i), err(:,i)] = RK34step(f, t(i-1), y(:,i-1), h(i));
err(:,i) = norm(err(:,i));

end

function dudt = lotka2(t,u)
a = 3;
b = 9;
c = 15;
d = 15;
dudt = [a*u(1)-b*u(1)*u(2); c*u(1)*u(2)-d*u(2)];
end










