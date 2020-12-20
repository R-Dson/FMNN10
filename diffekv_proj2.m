%% 1.1
ff = @(x,y) x.^5 - x.^3;
yprimprim = @(x, y) 20.*x.^(3) - 6.*x;
L = 2;
N = 1000; %interior points
alpha = 0; %då x = 0
beta = L^5-L.^3; %då x = L

x = linspace(L/(N+1), L-(L/(N+1)), N);
fvec = yprimprim(x);
y = twopBVP(fvec , alpha, beta, L, N);
y = [alpha y beta];
xfinal = linspace(0, L, N+2);

yreal = ff(xfinal);
error = zeros(size(y));
for n = 1:N+2
    error(n) = norm(yreal(n)-y(n));
end



%plot(xfinal, y);
hold on;
plot(xfinal, error);
xlabel('x');
ylabel('y');
%% 1.2 M
L = 10;
E = 1.9*10^11;
II = @(x, y) 10.^(-3).*(3-2.*cos(pi.*x./L).^12);
alpha = 0;
beta = 0;
N = 999;
q = @(x, y) -50*10^3.*x.^(0);
x = linspace(L/(N+1), L-(L/(N+1)), N);
xfinal = linspace(0, L, N+2);
fvec = q(x);

M = twopBVP(fvec , alpha, beta, L, N);

I = II(x);
fvec = M./(E.*I);

u = twopBVP(fvec, alpha, beta, L, N);
M = [alpha M beta];
u = [alpha u beta];

plot(xfinal, u);
xlabel('x');
ylabel('m');
u(501)

%% 2.1

Nend = 250;
Nstart = 20;
alleigs = zeros(3,Nend-Nstart+1);
for N = Nstart:Nend
    A = (diag(-2*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1));
    A(N,N) = A(N,N) + 4/3;
    A(N,N-1) = A(N,N-1) -1/3;
    A = (N+1)^2.*A;
    [eigf, eigs] = eig(A);
    eigs = eig(eigs);
    alleigs(1,N) = eigs(end);
    alleigs(2,N) = eigs(end-1);
    alleigs(3,N) = eigs(end-2);
end

for i = 1:Nstart-1
   alleigs(:,1) = [];
end
for i = 1:3
    alleigs(i,:) = alleigs(i,:) + (pi/2 + (i-1)*pi)^2;
end

x = linspace(Nstart, Nend, Nend-Nstart+1);
loglog(x, alleigs(1,:))
hold on;
loglog(x, alleigs(2,:))
loglog(x, alleigs(3,:))
loglog(x, 1./x.^2)
xlabel('N');
ylabel('error');
legend('first eigenvalue', 'second eigenvalue', 'third eigenvalue', 'reference plot');
grid on;


%% 2.1 igen
N = 499;
A = (diag(-2*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1));

A(N,N) = A(N,N) + 4/3;
A(N,N-1) = A(N,N-1) - 1/3;
%A(N-1,N-1)
%A(N+1,N) = 4/3 + 1;
%A(N+1, N-1) = -1/3 +1;
%A(N+1,N+1) = -2;

dx = 1/(N+1);
A = A./(dx^2);

[eigf, eigs] = eig(A);
eigs = diag(eigs);
[eigs, ind] = sort(eigs, 'descend'); %sorterar egenvärden i storleksordn
eigf = eigf(:,ind); %sorterar egenfunktioner i storleksordning

figure;
hold on;
x = linspace(0,1,N);
plot(x, eigf(:,1))
plot(x, eigf(:,2))
plot(x, eigf(:,3))

xlabel('x');
ylabel('y');
legend('first eigenmode', 'second eigenmode', 'third eigenmode');
grid on;
format long
eg = eigs(1)
eigs(2)
eigs(3)


%% 2.1 alt
exact = -2.467399062776626
N = 499;
M = N;
A = (diag(-2*ones(1,M)) + diag(1*ones(1,M-1),1) + diag(1*ones(1,M-1),-1));

A(M,M-1) = A(M,M-1) + 1;

dx = 1/(N);
A = A./(dx^2);

[eigf, eigs] = eig(A);
eigs = diag(eigs);
[eigs, ind] = sort(eigs, 'descend'); %sorterar egenvärden i storleksordn
eigf = eigf(:,ind); %sorterar egenfunktioner i storleksordning

format long
eigs(1)
eigs(2)
eigs(3)
diff = eigs(1) - exact

%% 2.2
%V = @(x,y) -2500*log((x-0.5).^2+1)
%V = @(x,y)  -9500*exp((x-0.5).^6) - 9500*(x-0.5).^2
%V = @(x,y) 1500*(0.3 - abs(x-0.5)); 
%V = @(x,y) 700*sin(pi.*x).^2;
%V = @(x,y) 3000*sin(4*pi.*x).^2;
V = @(x,y) (10^(-3))*(3-2.*cos(4*pi.*x).^12);

N = 399;
Vvec = V(linspace(1/(N+1),1-(1/(N+1)),N));


[eigv, eigf] = sturmliouville(Vvec, N);
x = linspace(0,1,N);

probdens = eigf.^2;


eigv
for i = 1:6
    probdens(:,i) = probdens(:,i)./(sum(probdens(:,i)));
    eigf(:,i) = eigf(:,i)./sum(probdens(:,i));
    probdens(:,i) = eigf(:,i).^2;
    eigf(:,i) = 500*eigf(:,i) - eigv(i);
end

figure
hold on

for i = 1:6
    probdens(:,i) = 20000*probdens(:,i) -eigv(i);
    plot(x, probdens(:,i))
end
xlabel('x')
ylabel('Energy level')

figure;
hold on;
for i = 1:6
    plot(x, eigf(:,i))
end

xlabel('x')
ylabel('Energy level')

%%

figure;
hold on
for i = 1:10
    probdens(:,i) = probdens(:,i) - eigv(i);
    plot(x, probdens(:,i))
end


function [eigval, eigfunc] = sturmliouville(Vvec, N)
    dx = 1/(N+1);
    A = (diag(-2*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1));
    A = A - dx^2.*diag(Vvec);
    
    A = (1./dx).^2.*A;
    [eigf, eigs] = eig(A);
    eigs = diag(eigs);
    [eigs, ind] = sort(eigs, 'descend'); %sorterar egenvärden i storleksordn
    eigf = eigf(:,ind); %sorterar egenfunktioner i storleksordning
    eigval = eigs(1:6);
    eigfunc = eigf(:,1:6);
end


function y = twopBVP(fvec, alpha, beta, L, N)
    dx = L/(N+1);
    A = diag(-2*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1);
    fvec(1) = fvec(1) -alpha/dx^2;
    fvec(end) = fvec(end) -beta/dx^2;
    y = dx^2.*fvec/A;
end