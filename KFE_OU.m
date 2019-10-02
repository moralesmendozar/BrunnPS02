function [p]=KFE_OU(ttheta,ssigma,xbar,x_grid,t_grid,p0,p_xmin,p_xmax,method,firstD)

% This function solves the Kolmogorov Forward Equation for a variable with
% a Ornstein-Uhlenbeck process of the form:
% dXt=ttheta(xbar-Xt).dt+ssigma.dZt
%
% Inputs:
% 
% 1. ttheta - "autoregressive parameter"
% 2. ssigma - variance of the shock
% 3. xbar - mean of the diffusion
% 4. x_grid - grid of x that we're going to use to numerically solve the
%    PDE (column vector)
% 5. t_grid - grid for time (row vector)
% 6. p0 - p(x,0) initial density on x (should have the same dimension as
%    x_grid and be a column vector)
% 7. p_xmin - p(xmin,.) boundary condition for the lowest x in the grid
%    (should have the same dimension as t_grid and be a row vector)
% 8. p_xmax - p(xmax,.) boundary condition for the largest x in the grid 
%    (should have the same dimension as t_grid and be a row vector)
% 9. method - =1 for Explicit Euler Method; =2 for Implicit Euler Method
% 10. firstD - first derivatives with respect to x: =1 for right
%     difference; =2 for central difference; =3 for left difference; =4 for
%     upwind scheme
%
% Outputs:
%
% p - p(x,t) a matrix which is nx-by-nt with the distribution of X over
% time

if t_grid(1)~=0
    fprintf('error: t_grid should start with t=0')
    return
end

nx=size(x_grid,1);
nt=size(t_grid,2);


if size(p0,1)~=nx
    fprintf('error: size of p0 is not equal to nx')
    return
elseif size(p_xmin,2)~=nt
    fprintf('error: size of p_xmin is not equal to nt')
    return
elseif size(p_xmax,2)~=nt
    fprintf('error: size of p_xmax is not equal to nt')
    return
end

p=zeros(nx,nt);
p(:,1)=p0;
p(1,:)=p_xmin;
p(nx,:)=p_xmax;
    
    
if method==1
    

    
    for t=1:nt-1
        for x=2:nx-1
            if firstD==1
                dpdx(x,t)=(p(x+1,t)-p(x,t))/(x_grid(x+1)-x_grid(x));
            elseif firstD==2
                dpdx(x,t)=(p(x+1,t)-p(x-1,t))/(x_grid(x+1)-x_grid(x-1));
            elseif firstD==3
                dpdx(x,t)=(p(x,t)-p(x-1,t))/(x_grid(x)-x_grid(x-1));
            elseif firstD==4
                if x_grid(x)>xbar
                    dpdx(x,t)=(p(x+1,t)-p(x,t))/(x_grid(x+1)-x_grid(x));
                else
                    dpdx(x,t)=(p(x,t)-p(x-1,t))/(x_grid(x)-x_grid(x-1));
                end
            else
                fprintf('error: inappropriate input for first derivative method');
                return
            end
            d2pdx(x,t)=(p(x+1,t)+p(x-1,t)-2*p(x,t))/(x_grid(x+1)-x_grid(x))^2;
            p(x,t+1)=(1+ttheta*(t_grid(t+1)-t_grid(t)))*p(x,t)+(t_grid(t+1)-t_grid(t))*...
            (ttheta*(x_grid(x)-xbar)*dpdx(x,t)+((ssigma^2)/2)*d2pdx(x,t));
        end
    end
       
    

elseif method==2
        
    
    for t=1:nt-1
        b=zeros(nx,1);
        b(1,1)=p(1,t+1);
        b(nx,1)=p(nx,t+1);
        b(2:nx-1,1)=p(2:nx-1,t);
        A=spalloc(nx,nx,2+3*(nx-2));
        A(1,1)=1;
        A(nx,nx)=1;
        for x=2:nx-1
            ind=[x-1,x,x+1];
            if firstD==1
                aux=[-((t_grid(t+1)-t_grid(t))*ssigma^2)/(2*(x_grid(x+1)-x_grid(x))^2),...
                1-(t_grid(t+1)-t_grid(t))*ttheta+((t_grid(t+1)-t_grid(t))*ttheta*(x_grid(x)-xbar))/...
                (x_grid(x+1)-x_grid(x))+((t_grid(t+1)-t_grid(t))*ssigma^2)/(x_grid(x+1)-x_grid(x))^2,...
                -((t_grid(t+1)-t_grid(t))*ttheta*(x_grid(x)-xbar))/(x_grid(x+1)-x_grid(x))-...
                ((t_grid(t+1)-t_grid(t))*ssigma^2)/(2*(x_grid(x+1)-x_grid(x))^2)];
            elseif firstD==2
                aux=[-((t_grid(t+1)-t_grid(t))*ssigma^2)/(2*(x_grid(x+1)-x_grid(x))^2)+...
                ((t_grid(t+1)-t_grid(t))*ttheta*(x_grid(x)-xbar))/(x_grid(x+1)-x_grid(x-1)),...
                1-(t_grid(t+1)-t_grid(t))*ttheta+((t_grid(t+1)-t_grid(t))*ssigma^2)/...
                (x_grid(x+1)-x_grid(x))^2,-((t_grid(t+1)-t_grid(t))*ttheta*(x_grid(x)-xbar))/...
                (x_grid(x+1)-x_grid(x-1))-((t_grid(t+1)-t_grid(t))*ssigma^2)/...
                (2*(x_grid(x+1)-x_grid(x))^2)];
            elseif firstD==3
                aux=[-((t_grid(t+1)-t_grid(t))*ssigma^2)/(2*(x_grid(x+1)-x_grid(x))^2)+...
                ((t_grid(t+1)-t_grid(t))*ttheta*(x_grid(x)-xbar))/(x_grid(x)-x_grid(x-1)),...
                1-(t_grid(t+1)-t_grid(t))*ttheta+((t_grid(t+1)-t_grid(t))*ssigma^2)/...
                (x_grid(x+1)-x_grid(x))^2-((t_grid(t+1)-t_grid(t))*ttheta*(x_grid(x)-xbar))/...
                (x_grid(x)-x_grid(x-1)),-((t_grid(t+1)-t_grid(t))*ssigma^2)/...
                (2*(x_grid(x+1)-x_grid(x))^2)];
            else
                fprintf('error: inappropriate input for first derivative method')
                return
            end
               A(x,:)=sparse(1,ind,aux,1,nx,3);
        end
        p(:,t+1)=A\b;
    end

else
    fprintf('error: method input should equal 1 or 2')
    return

end


end

