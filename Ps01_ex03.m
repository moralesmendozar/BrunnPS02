% UPENN
% Brunnermeier online Course / Princeton
% September 21, 2019.
% Problem set 01, ex 03
%close all;
clear;
clc;
n = 2000;
options = optimset('Display', 'off');
xx = linspace(0,10,n);
y_explicitEuler = ones(n,1);
y_implicitEuler = ones(n,1);
yreal = exp(xx);
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