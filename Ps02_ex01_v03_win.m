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
xstep = 0.01; %0.1  %then 0.05
x_grid = xmin:xstep:xmax;
n = length(x_grid);
tmin = 0;
tmax = 1;
tstep = 0.0001;%0.1%0.005;
t_grid = tmin:tstep:tmax;
T = length(t_grid);
xx = x_grid'*ones(1,T);
tt = ones(n,1)*t_grid;
% parameters:
theta = 0;
sigma = 0.33;
xbar = 0;
m0 = -3;
v0 = 0.33;
p0 = normpdf(x_grid,m0,v0);
pN = p0;
if theta == 0
    v =@(t) v0;
    m =@(t) m0;
else
    v =@(t) v0*exp(-2*theta.*t) + (1-exp(-2*theta.*t)).*sigma^2./(2*theta);
    m =@(t) m0*exp(-theta.*t) + (1-exp(-theta.*t)).*xbar;
end
%Analytical solution to the PDE
p_real =@(x,t) normpdf( (x-m(t))./(sqrt(v(t))))./(sqrt(v(t)));
preal = p_real(xx,tt);
%Explicit, ImplicitEuler Methods with forward, backward, etc.
display('Running explicitEuler with central Derivative...')
display('Running time = ...')
tic
p_explicitEuler = fn_KolmogForwEqn(xx,tt,options, theta, sigma, xbar, p0, pN, 1,2);
toc
display('Running implicitEuler with central Derivative...')
display('Running time = ...')
tic
p_implicitEuler = fn_KolmogForwEqn(xx,tt,options, theta, sigma, xbar, p0, pN, 2,2);
toc
tic
%p_matlabSolver = fn_KolmogForwEqn(xx,tt,options, theta, sigma, xbar, p0, pN, 2,1);
display('Running matlabSolver with central Derivative...')
display('Running time = ...')
%3D plot p_explicit vs preal

columna = 10;
figure(1)
plot(x_grid,p_explicitEuler(:,columna),'--r','Linewidth',2)
xlabel('x')
texty = ['P(x,t=',num2str(columna),')'];
ylabel(texty)
hold on;
plot(x_grid,p_implicitEuler(:,columna),'-*b','Linewidth',0.01)
%plot(x_grid,p_matlabSolver(:,columna),'--g','Linewidth',2)
plot(x_grid,preal(:,columna),'-k')
%legend('Explicit Euler','Implicit Euler','Matlab Solution','Analytical solution');
legend('Explicit Euler','Implicit Euler','Analytical solution')
hold off;





