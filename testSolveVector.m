options = optimset('Display', 'off');
clear
clc;
n = 10000;
A = round(rand(n)*10);
A = A*A';
ee = eig(A);
pert = min(ee);
%if (pert<=0)
if (det(A)<=0)
    A = A + pert*eye(n);
end
b = rand(n,1)*100;

tic
x = A\b;
toc
sum(abs(A*x - b))
%y =@(x) A*x - b;
%init = [0,0];
%[ysolve,fval] = fsolve(f,init,options);
%ysolve
%y(ysolve)
%%
AA = A;
for ii =1:n
    for jj = 1:n
        if ii>=jj+2 || ii <= jj-2
            AA(ii,jj) = 0;
        end
    end
end
tic
x = AA\b;
toc
sum(abs(AA*x - b))