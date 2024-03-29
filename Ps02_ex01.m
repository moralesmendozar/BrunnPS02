% UPENN
% Brunnermeier online Course / Princeton
% September 21, 2019.
% Problem set 01, ex 03
%close all;
clear;
clc;

%n = 100;
%T = 200;
options = optimset('Display', 'off');
xmin = -5;
xmax = 5;
xstep = 0.1; %0.1
x_grid = xmin:xstep:xmax;
n = length(x_grid);
tmin = 0;
tmax = 2;
tstep = 0.1;%0.005;
t_grid = tmin:tstep:tmax;
T = length(t_grid);
xx = x_grid'*ones(1,T);
tt = ones(n,1)*t_grid;
% parameters:
theta = 0;
sigma = 1;
xbar = 0;
m0 = 0;
v0 = 0.1;
v =@(t) v0*exp(-2*theta.*t) + (1-exp(-2*theta.*t)).*sigma^2./(2*theta);
m =@(t) m0*exp(-theta.*t) + (1-exp(-theta.*t)).*xbar;
%Analytical solution to the PDE
p_real =@(x,t) normpdf( (x-m(t))./(sqrt(v(t))))./(sqrt(v(t)));
preal = p_real(xx,tt);
%Explicit, ImplicitEuler Methods with forward, backward, etc.
p_explicitEuler = fn_KolmogForwEqn(xx,tt,options, theta, sigma, xbar, 1,1);
p_implicitEuler = fn_KolmogForwEqn(xx,tt,options, theta, sigma, xbar, 2,1);
%3D plot p_explicit vs preal










g=@(x,y) y;   %this is y' = y   %for exercise a)

for counter = 2:n
    xi = xx(counter-1);
    yi  = y_explicitEuler(counter-1);
    y_explicitEuler(counter) = yi + g(xi,yi)*(xx(counter)-xi);
    y_implicitEuler(counter) = yi/(1-(xx(counter) - xi));
end

%use matlab to solve
syms z(t);
ode = diff(z,t) == z;
cond = z(0) == 1;
ySol(t) = dsolve(ode,cond);


y_explicitEuler = 


figure(1)
plot(xx,y_explicitEuler,'--r','Linewidth',3)
xlabel('x')
ylabel('y')
hold on;
plot(xx,y_implicitEuler,'-*b','Linewidth',0.01)
plot(xx,ySol(xx),'--g','Linewidth',2)
plot(xx,yreal,'-k')
legend('Explicit Euler','Implicit Euler','Matlab Solution','Analytical solution');
hold off;

%% ex 03_b
xx = linspace(0,10,n);
y_explicitEuler = ones(n,1);
y_implicitEuler = ones(n,1);
yreal = 1./(1-sin(xx.^2)./2);
g=@(x,y) x*cos(x^2)*y^2;   %this is y'
for counter = 2:n
    xi = xx(counter-1);
    yi  = y_explicitEuler(counter-1);
    y_explicitEuler(counter) = yi + g(xi,yi)*(xx(counter)-xi);
    f=@(y) (y-yi)/(xx(counter)-xi) - g(xx(counter),y);
    [ysolve,fval] = fsolve(@(y)f(y),yi,options);
    y_implicitEuler(counter) = ysolve;
end
%use matlab to solve
syms z(t);
ode = diff(z,t) == t*cos(t^2)*z^2;
cond = z(0) == 1;
ySol(t) = dsolve(ode,cond);

figure(2)
plot(xx,y_explicitEuler,'--r','Linewidth',3)
xlabel('x')
ylabel('y')
hold on;
plot(xx,y_implicitEuler,'-*b')
plot(xx,ySol(xx),'--g','Linewidth',2)
plot(xx,yreal,'-k')
legend('Explicit Euler','Implicit Euler','Matlab Solution','Analytical solution');
hold off;

%% ex 03_c
xx = linspace(0,10,n);
y_explicitEuler = ones(n,1);
y_implicitEuler = ones(n,1);
v = zeros(n,1);  % dy/dx = v,   dv/dt = d^2y/dx^2 = -y
yreal = cos(xx);
g=@(x,y) -y;   %this is y'' = -y

for counter = 2:n
    xi = xx(counter-1);
    yi  = y_explicitEuler(counter-1);
    vi  = v(counter-1);
    v(counter) = vi - yi*(xx(counter)-xi);
    y_explicitEuler(counter) = yi + v(counter)*(xx(counter)-xi);
    %f=@(y) (y-yi)/(xx(counter)-xi) - g(xx(counter),y);
    %[ysolve,fval] = fsolve(@(y)f(y),yi,options);
    %y_implicitEuler(counter) = ysolve;
    y_implicitEuler(counter) = y_explicitEuler(counter);
end
%use matlab to solve
syms z(t);
Dz = diff(z);
ode = diff(z,t,2) == -z;
cond1 = z(0) == 1;
cond2 = Dz(0) == 0;
conds = [cond1 cond2];
ySol(t) = dsolve(ode,conds);

figure(3)
plot(xx,y_explicitEuler,'--r','Linewidth',3)
xlabel('x')
ylabel('y')
hold on;
%plot(xx,y_implicitEuler,'-*b')
plot(xx,ySol(xx),'--g','Linewidth',2)
plot(xx,yreal,'-k')
legend('Explicit Euler','Matlab Solution','Analytical solution');
hold off;