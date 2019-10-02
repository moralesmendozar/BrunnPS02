% UPENN
% Brunnermeier online Course / Princeton
% September 21, 2019.
% Problem set 01, ex 03
%close all;

function [P] = fn_KolmogForwEqn(xx,tt,options, theta, sigma, xbar, p0, pN, type,derivtype)
% This function solves p(x,t) for the diff eqn
%   d p(x,t)/dt = -theta*(x - xbar)*dp/dx + theta*p +sigma^2*d2p/dx2
%  inputs are:
%       xx, tt vectors containing xs and tts, in the form
%           xx = x_grid'*ones(1,T)   tt = ones(n,1)*t_grid
%       options for optim, 
%       theta, sigma, xbar parameteres of the model
%       p0 initial distribution (e.g. a normal one)
%  type == 1 = explicit euler, 2 = implicit Euler, 3 = Matlab method
%  derivtype == 1, central, 2= left, 3=right, 4= upwind
    [n T] = size(xx);
    t_grid = tt(1,:);
    x_grid = xx(:,1);
    P = ones(n,T); %will save the solution
    P(:,1) = p0;
    P(:,T) = pN;
    P(1,:) = P(1,1);
    P(n,:) = P(n,T);
    if type ==1
        %Explicit Euler Method
        if derivtype == 1 %center, (Explicit)
            %center, (Explicit)
            for t = 2:T
                deltat = t_grid(t) - t_grid(t-1);
                for ii = 2:(n-1)
                    deltaxd = x_grid(ii+1) - x_grid(ii-1);
                    deltax =  x_grid(ii) - x_grid(ii-1);
                    pstatic = theta *P(ii,t-1);
                    pde_firstDeriv = theta*( x_grid(ii) - xbar)*(P(ii+1,t-1)-P(ii-1,t-1))/(deltaxd);
                    pde_secondDeriv = sigma^2/2*(P(ii+1,t-1)+P(ii-1,t-1)-2*P(ii,t-1))/deltax^2;
                    pde_it = pstatic + pde_firstDeriv + pde_secondDeriv;
                    P(ii,t) = P(ii,t-1) + deltat*(pde_it);
                end
            end
        elseif derivtype == 2 %left, (Explicit)
            %left, (Explicit)
            for t = 2:T
                deltat = t_grid(t) - t_grid(t-1);
                for ii = 2:(n-1)
                    deltaxd = x_grid(ii) - x_grid(ii-1);
                    deltax =  x_grid(ii) - x_grid(ii-1);
                    pstatic = theta *P(ii,t-1);
                    pde_firstDeriv = theta*( x_grid(ii) - xbar)*(P(ii,t-1)-P(ii-1,t-1))/(deltaxd);
                    pde_secondDeriv = sigma^2/2*(P(ii+1,t-1)+P(ii-1,t-1)-2*P(ii,t-1))/(deltax)^2;
                    pde_it = pstatic + pde_firstDeriv + pde_secondDeriv;
                    P(ii,t) = P(ii,t-1) + deltat*(pde_it);
                end
            end
        elseif derivtype == 3%right, (Explicit)
            %right, (Explicit)
            for t = 2:T
                deltat = t_grid(t) - t_grid(t-1);
                for ii = 2:(n-1)
                    deltaxd = x_grid(ii+1) - x_grid(ii);
                    deltax =  x_grid(ii) - x_grid(ii-1);
                    pstatic = theta *P(ii,t-1);
                    pde_firstDeriv = theta*( x_grid(ii) - xbar)*(P(ii+1,t-1)-P(ii,t-1))/(deltaxd);
                    pde_secondDeriv = sigma^2/2*(P(ii+1,t-1)+P(ii-1,t-1)-2*P(ii,t-1))/(deltax)^2;
                    pde_it = pstatic + pde_firstDeriv + pde_secondDeriv;
                    P(ii,t) = P(ii,t-1) + deltat*(pde_it);
                end
            end
        elseif derivtype == 4%upwind, (Explicit)
            %upwind, (Explicit)
            for t = 2:T
                deltat = t_grid(t) - t_grid(t-1);
                for ii = 2:(n-1)
                    pstatic = theta *P(ii,t-1);
                    coeffOfDx = theta*( x_grid(ii) - xbar);
                    if coeffOfDx >= 0 
                        deltax = x_grid(ii+1) - x_grid(ii);
                        derivRight = (P(ii+1,t-1)-P(ii,t-1))/(deltax);
                        pde_firstDeriv = coeffOfDx*derivRight;
                    else
                        deltax = x_grid(ii) - x_grid(ii-1);
                        derivLeft = (P(ii,t-1)-P(ii-1,t-1))/(deltax);
                        pde_firstDeriv = coeffOfDx*derivLeft;
                    end
                    deltax =  x_grid(ii) - x_grid(ii-1);
                    pde_secondDeriv = sigma^2/2*(P(ii+1,t-1)+P(ii-1,t-1)-2*P(ii,t-1))/(deltax)^2;
                    pde_it = pstatic + pde_firstDeriv + pde_secondDeriv;
                    P(ii,t) = P(ii,t-1) + deltat*(pde_it);
                end
            end
        else
            display('Error, derivType not found')
        end
    elseif type ==2
        %Implicit Euler Method
        if derivtype == 1 %center, implicit
            %center, implicit
            for t = 2:T
                deltat = t_grid(t) - t_grid(t-1);
                A = zeros(n-2);
                b = zeros(n-2,1);
                deltax1d = x_grid(3)-x_grid(1);
                deltax1 = x_grid(2) - x_grid(1);
                A(1,1) = (1+ deltat*sigma^2/deltax1^2 -deltat*theta);
                A(1,2) = -deltat*theta*(xx(2,t)-xbar)/deltax1d - deltat*sigma^2/(2*deltax1^2);
                deltaxnd = x_grid(n)-x_grid(n-2);
                deltaxn = x_grid(n-1)-x_grid(n-2);
                A(n-2,n-3) = deltat*theta*(xx(n-1,t)-xbar)/deltaxnd - deltat*sigma^2/(2*deltaxn^2);
                A(n-2,n-2) = (1+ deltat*sigma^2/deltaxn^2 -deltat*theta);
                b(1) = P(2,t-1) + P(1,t)*(-deltat*theta*(xx(2,t)-xbar)/deltax1d+deltat*sigma^2/(2*deltax1^2));
                b(n-2) = P(n-1,t-1) + P(n,t)*(deltat*sigma^2/(2*deltaxn^2) + deltat*theta*(xx(n-1,t)-xbar)/deltaxnd);
                for ii = 2:(n-3)
                    deltaxid = x_grid(ii+2) - x_grid(ii);
                    deltaxi = x_grid(ii+1) - x_grid(ii);
                    A(ii,ii-1) = deltat*theta*(xx(ii+1,t)-xbar)/deltaxid - deltat*sigma^2/(2*deltaxi^2);
                    A(ii,ii) = (1+ deltat*sigma^2/deltaxi^2 -deltat*theta);
                    A(ii,ii+1) = -deltat*theta*(xx(ii+1,t)-xbar)/deltaxid - deltat*sigma^2/(2*deltaxi^2);
                    b(ii) = P(ii+1,t-1);
                end
                Psolvetemp = A\b;
                P(2:(n-1),t) = Psolvetemp;
            end
        elseif derivtype == 2 %left, implicit
            %left, Implicit
            for t = 2:T
                deltat = t_grid(t) - t_grid(t-1);
                A = zeros(n-2);
                b = zeros(n-2,1);
                deltax1d = x_grid(2)-x_grid(1);
                deltax1 = x_grid(2)-x_grid(1);
                A(1,1) = (1-deltat*theta*(xx(2,t)-xbar)/deltax1d + deltat*sigma^2/deltax1^2 -deltat*theta);
                A(1,2) = - deltat*sigma^2/(2*deltax1^2);
                deltaxnd = x_grid(n-1)-x_grid(n-2);
                deltaxn = x_grid(n-1)-x_grid(n-2);
                A(n-2,n-3) = deltat*theta*(xx(n-1,t)-xbar)/deltaxnd - deltat*sigma^2/(2*deltaxn^2);
                A(n-2,n-2) = (1-deltat*theta*(xx(n-1,t)-xbar)/deltaxnd + deltat*sigma^2/deltaxn^2 -deltat*theta);
                b(1) = P(2,t-1) + P(1,t)*(-deltat*theta*(xx(2,t)-xbar)/deltax1d+deltat*sigma^2/(2*deltax1^2));
                b(n-2) = P(n-1,t-1) + P(n,t)*(deltat*sigma^2/(2*deltaxn^2));
                for ii = 2:(n-3)
                    deltaxid = x_grid(ii+1) - x_grid(ii);
                    deltaxi = x_grid(ii+1) - x_grid(ii);
                    A(ii,ii-1) = deltat*theta*(xx(ii+1,t)-xbar)/deltaxid - deltat*sigma^2/(2*deltaxi^2);
                    A(ii,ii) = (1-deltat*theta*(xx(ii+1,t)-xbar)/deltaxid + deltat*sigma^2/deltaxi^2 -deltat*theta);
                    A(ii,ii+1) = - deltat*sigma^2/(2*deltaxi^2);
                    b(ii) = P(ii+1,t-1);
                end
                Psolvetemp = A\b;
                P(2:(n-1),t) = Psolvetemp;
            end
        elseif derivtype == 3%right, implicit
            %right, implicit
            for t = 2:T
                deltat = t_grid(t) - t_grid(t-1);
                A = zeros(n-2);
                b = zeros(n-2,1);
                deltax1d = x_grid(3)-x_grid(2);
                deltax1 = x_grid(2)-x_grid(1);
                A(1,1) = (1 +deltat*theta*(xx(2,t)-xbar)/deltax1d + deltat*sigma^2/deltax1^2 -deltat*theta);
                A(1,2) = - deltat*sigma^2/(2*deltax1^2)-deltat*theta*(xx(2,t)-xbar)/deltax1d;
                deltaxnd = x_grid(n)-x_grid(n-1);
                deltaxn = x_grid(n-1)-x_grid(n-2);
                A(n-2,n-3) = - deltat*sigma^2/(2*deltaxn^2);
                A(n-2,n-2) = (1 +deltat*theta*(xx(n-1,t)-xbar)/deltaxnd + deltat*sigma^2/deltaxn^2 -deltat*theta);
                b(1) = P(2,t-1) + P(1,t)*(deltat*sigma^2/(2*deltax1^2));
                b(n-2) = P(n-1,t-1) + P(n,t)*(deltat*sigma^2/(2*deltaxn^2)+deltat*theta*(xx(n-1,t)-xbar)/deltaxnd);
                for ii = 2:(n-3)
                    deltaxid = x_grid(ii+2) - x_grid(ii+1);
                    deltaxi = x_grid(ii+1) - x_grid(ii+0);
                    A(ii,ii-1) = - deltat*sigma^2/(2*deltaxi^2);
                    A(ii,ii) = (1 + deltat*theta*(xx(ii+1,t)-xbar)/deltaxid + deltat*sigma^2/deltaxi^2 -deltat*theta);
                    A(ii,ii+1) = - deltat*sigma^2/(2*deltaxi^2)-deltat*theta*(xx(ii+1,t)-xbar)/deltaxid;
                    b(ii) = P(ii+1,t-1);
                end
                Psolvetemp = A\b;
                P(2:(n-1),t) = Psolvetemp;
            end
        elseif derivtype == 4%upwind, (Implicit)
            %upwind, (Implicit)
            for t = 2:T
                deltat = t_grid(t) - t_grid(t-1);
                A = zeros(n-2);
                b = zeros(n-2,1);
                deltax1 = x_grid(2)-x_grid(1);
                coeffOfDx = deltat*theta*(xx(2,t)-xbar);
                if coeffOfDx >= 0
                    %Use right
                    deltax1d = x_grid(3)-x_grid(2);
                    A(1,1) = (1 +coeffOfDx/deltax1d + deltat*sigma^2/deltax1^2 -deltat*theta);
                    A(1,2) = - deltat*sigma^2/(2*deltax1^2)-coeffOfDx/deltax1d;
                    b(1) = P(2,t-1) + P(1,t)*(deltat*sigma^2/(2*deltax1^2));
                else
                    %use left
                    deltax1d = x_grid(2)-x_grid(1);
                    A(1,1) = (1-coeffOfDx/deltax1d + deltat*sigma^2/deltax1^2 -deltat*theta);
                    A(1,2) = - deltat*sigma^2/(2*deltax1^2);
                    b(1) = P(2,t-1) + P(1,t)*(-coeffOfDx/deltax1d+deltat*sigma^2/(2*deltax1^2));
                end
                deltaxn = x_grid(n-1)-x_grid(n-2);
                coeffOfDx = deltat*theta*(xx(n-1,t)-xbar);
                if coeffOfDx >= 0
                    deltaxnd = x_grid(n)-x_grid(n-1);
                    %Use right
                    A(n-2,n-3) = - deltat*sigma^2/(2*deltaxn^2);
                    A(n-2,n-2) = (1 +coeffOfDx/deltaxnd + deltat*sigma^2/deltaxn^2 -deltat*theta);
                    b(n-2) = P(n-1,t-1) + P(n,t)*(deltat*sigma^2/(2*deltaxn^2)+coeffOfDx/deltaxnd);
                else
                    %use left
                    deltaxnd = x_grid(n-1)-x_grid(n-2);
                    A(n-2,n-3) = coeffOfDx/deltaxnd - deltat*sigma^2/(2*deltaxn^2);
                    A(n-2,n-2) = (1-coeffOfDx/deltaxnd + deltat*sigma^2/deltaxn^2 -deltat*theta);
                    b(n-2) = P(n-1,t-1) + P(n,t)*(deltat*sigma^2/(2*deltaxn^2));
                end
                for ii = 2:(n-3)
                    coeffOfDx = deltat*theta*(xx(ii+1,t)-xbar);
                    if coeffOfDx >= 0
                        %Use right
                        deltaxid = x_grid(ii+2) - x_grid(ii+1);
                        deltaxi = x_grid(ii+1) - x_grid(ii+0);
                        A(ii,ii-1) = - deltat*sigma^2/(2*deltaxi^2);
                        A(ii,ii) = (1 + coeffOfDx/deltaxid + deltat*sigma^2/deltaxi^2 -deltat*theta);
                        A(ii,ii+1) = - deltat*sigma^2/(2*deltaxi^2)-coeffOfDx/deltaxid;
                        b(ii) = P(ii+1,t-1);
                    else
                        %use left
                        deltaxid = x_grid(ii+1) - x_grid(ii);
                        deltaxi = x_grid(ii+1) - x_grid(ii+0);
                        A(ii,ii-1) = coeffOfDx/deltaxid - deltat*sigma^2/(2*deltaxi^2);
                        A(ii,ii) = (1-coeffOfDx/deltaxid + deltat*sigma^2/deltaxi^2 -deltat*theta);
                        A(ii,ii+1) = - deltat*sigma^2/(2*deltaxi^2);
                        b(ii) = P(ii+1,t-1);
                    end
                end
                Psolvetemp = A\b;
                P(2:(n-1),t) = Psolvetemp;
            end
        else
            display('Error, derivType not found')
        end
    elseif type ==3
        %use matlab to solve
        syms z(t);
        ode = diff(z,t) == t*cos(t^2)*z^2;
        cond = z(0) == 1;
        ySol(t) = dsolve(ode,cond);
    else
        display('Error, type Method not found')
    end
    %the function returns P in the end.    
end