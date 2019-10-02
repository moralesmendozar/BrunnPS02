% UPENN
% Brunnermeier online Course / Princeton
% September 21-okt1, 2019.
% Problem set 02, ex 03
close all;
clear;
clc;
exercise = '04';
options = optimset('Display', 'off');
xmin = -5;
xmax = 5;
% x_grid = xmin;
% xstep = rand()/10;
% if xstep < 0.001
%     xstep = 0.01;
% end
% xnext = xmin + xstep;
% while xnext < xmax
%     x_grid = [x_grid xnext];
%     %xstep = rand(); %0.1  %then 0.05
%     xstep = rand()/10;
%     if xstep < 0.001
%         xstep = 0.01;
%     end
%     xnext = xnext + xstep;
% end
% x_grid = [x_grid xmax];
xstep = 0.1; %0.1, 0.05   
x_grid = xmin:xstep:xmax;
n = length(x_grid);
tmin = 0;
tmax = 1;  %2,  1
tstep = 0.005;%0.1%0.005;
t_grid = tmin:tstep:tmax;
T = length(t_grid);
xx = x_grid'*ones(1,T);
tt = ones(n,1)*t_grid;
% parameters:
theta =3;  %0 = 0. 001,   3, 0.5
sigma = 0.33; %1, 0.33,   0
xbar = 0;
m0 = -3;  %0, -3,  0
v0 = 0.1; %0.1  0.33,  0.1
p0 = normpdf(x_grid,m0,sqrt(v0));
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
for derivtype = 1:4
    %Explicit, ImplicitEuler Methods with forward, backward, etc.
    if derivtype == 1
        word = 'Central';
    elseif derivtype ==2
        word = 'Left';
    elseif derivtype ==3
        word = 'Right';
    elseif derivtype ==4
        word = 'Upwind';
    else
        word = 'no such method';
        display(word)
    end
    display(['Running explicitEuler with ',word, ' Derivative...'])
    %display('Running time = ...')
    tic
    p_explicitEuler = fn_KolmogForwEqn(xx,tt,options, theta, sigma, xbar, p0, pN, 1,derivtype);
    toc
    display(['Running implicitEuler with ',word, ' Derivative...'])
    %display('Running time = ...')
    tic
    p_implicitEuler = fn_KolmogForwEqn(xx,tt,options, theta, sigma, xbar, p0, pN, 2,derivtype);
    toc
    %tic
    %display(['Running matlabSolver with ',word, ' Derivative...'])
    %p_matlabSolver = fn_KolmogForwEqn(xx,tt,options, theta, sigma, xbar, p0, pN, 3,1);
    % toc
%end

% Plot of P(x,t=columna) for explicit, implicit, et.al.
    counterFig = 0;
%for derivtype = 1:4
%     close all;
%     if derivtype == 1
%         word = 'Central';
%     elseif derivtype ==2
%         word = 'Left';
%     elseif derivtype ==3
%         word = 'Right';
%     elseif derivtype ==4
%         word = 'Upwind';
%     else
%         word = 'no such method';
%         display(word)
%     end
    counterFig = counterFig + 1;
    columna = 15;
    %2D plot p_explicit vs preal
    countf2 = counterFig;
    figure(countf2)
    plot(x_grid,p_explicitEuler(:,columna),'--r','Linewidth',2)
    xlabel('x')
    texty = ['P(x,t=',num2str(columna),')'];
    ylabel(texty)
    hold on;
    plot(x_grid,p_implicitEuler(:,columna),'-*b','Linewidth',0.01)
    plot(x_grid,preal(:,columna),'-k')
    legend('Explicit Euler','Implicit Euler','Analytical solution')
    title([word,' Derivative Method'])
    hold off;
    plotName = ['plots\ex',exercise,'_P(x,t=15)',num2str(derivtype)];
    print(plotName,'-dpng')
%end

%% 3D plot log_error of  |p_explicit - preal|
%counterFig = 0;
%for derivtype = 1:4
    close all;
%     if derivtype == 1
%         word = 'Central';
%     elseif derivtype ==2
%         word = 'Left';
%     elseif derivtype ==3
%         word = 'Right';
%     elseif derivtype ==4
%         word = 'Upwind';
%     else
%         word = 'no such method';
%         display(word)
%     end
    counterFig = counterFig + 1;
    %3D plot log_error of  p_explicit vs preal
    [X,Y] = meshgrid(t_grid,x_grid);
    figure(counterFig)
    surf(X,Y,log(abs(p_explicitEuler-preal)) )
    xlabel('t')
    ylabel('x')
    textz = ['log( |P(x,t) - Preal(x,t)| ) '];
    zlabel(textz)
    title(['3D plot for Explicit ', word,' Derivative Method'])
    plotName = ['plots\ex',exercise,'_3dPlot_AbsError_Explicit',num2str(derivtype)];
    print(plotName,'-dpng')
    %3D plot log_error of  p_implicitEuler vs preal
    counterFig = counterFig + 1;
    figure(counterFig)
    surf(X,Y,log(abs(p_implicitEuler-preal)) )
    xlabel('t')
    ylabel('x')
    textz = ['log( |P(x,t) - Preal(x,t)| ) '];
    zlabel(textz)
    title(['3D plot for Implicit ', word,' Derivative Method'])
    plotName = ['plots\ex',exercise,'_3dPlot_AbsError_Implicit',num2str(derivtype)];
    print(plotName,'-dpng')
%end

%% 3D plot log_error of  |p_explicit - preal| / | preal|
%counterFig = 0;
%for derivtype = 1:4
%   close all;
%     if derivtype == 1
%         word = 'Central';
%     elseif derivtype ==2
%         word = 'Left';
%     elseif derivtype ==3
%         word = 'Right';
%     elseif derivtype ==4
%         word = 'Upwind';
%     else
%         word = 'no such method';
%         display(word)
%     end
    counterFig = counterFig + 1;
    %3D plot log_error of relative p_explicit vs preal
    [X,Y] = meshgrid(t_grid,x_grid);
    figure(counterFig)
    Z_expl_errors = log(abs(p_explicitEuler-preal)./ abs( preal));
    surf(X,Y,Z_expl_errors)
    xlabel('t')
    ylabel('x')
    textz = ['log( |P(x,t) - Preal(x,t)| / |preal| ) '];
    zlabel(textz)
    title(['3D plot for Explicit ', word,' Derivative Method'])
    plotName = ['plots\ex',exercise,'_3dPlot_RelError_Explicit',num2str(derivtype)];
    print(plotName,'-dpng')
    
    %3D plot log_error of relative p_implicit vs preal
    counterFig = counterFig + 1;
    figure(counterFig)
    Z_impl_errors = log(abs(p_implicitEuler-preal)./ abs( preal));
    surf(X,Y,Z_impl_errors)
    xlabel('t')
    ylabel('x')
    textz = ['log( |P(x,t) - Preal(x,t)| / |preal| ) '];
    zlabel(textz)
    title(['3D plot for Implicit ', word,' Derivative Method'])
    plotName = ['plots\ex',exercise,'_3dPlot_RelError_Implicit',num2str(derivtype)];
    print(plotName,'-dpng')
    
    %Plot the maximum relative Error along the space dimension
    countf2 = counterFig;
    figure(20+ countf2)
    plot(t_grid,max(Z_expl_errors),'-r','Linewidth',1)
    xlabel('t')
    texty = [' max( log( |P(x,t) - Preal(x,t)| / |preal| )'];
    ylabel(texty)
    hold on;
    plot(t_grid,max(Z_impl_errors),'-k','Linewidth',1)
    title(['Max relative error along x axis w/ ', word,' Derivative Method'])
    plotName = ['plots\ex',exercise,'_max(abs_error)',num2str(derivtype)];
    hold off;
    legend({'Explicit Euler','Implicit Euler'},'Location','Southeast')
    print(plotName,'-dpng')
    %Plot the maximum relative Error along t for x
    countf2 = counterFig;
    figure(20+ countf2)
    plot(x_grid,max(Z_expl_errors'),'-r','Linewidth',1)
    xlabel('x')
    texty = [' max( log( |P(x,t) - Preal(x,t)| / |preal| )'];
    ylabel(texty)
    hold on;
    plot(x_grid,max(Z_impl_errors'),'-k','Linewidth',1)
    title(['Max relative error along t axis w/ ', word,' Derivative Method'])
    plotName = ['plots\ex',exercise,'_max(abs_error)_x',num2str(derivtype)];
    hold off;
    legend({'Explicit Euler','Implicit Euler'},'Location','Southeast')
    print(plotName,'-dpng')
    close all;
end


