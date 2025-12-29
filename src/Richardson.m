function [x,res] = Richardson(A,b,Tol)
%% Risolutore di sistemi lineari per mezzo del metodo di Richardson
x = zeros(size(b));
alpha = 2/(eigs(A,1)+eigs(A,1,'smallestabs'));     %% alpha ottimale, 2 diviso la somma tra lambda massimo e minimo
maxiter = 1000;
i=0;
res = zeros(1,maxiter);
r = norm(A*x-b);                                   %% residuo iniziale con x=b
res(1) = r;
while r>Tol && i<maxiter                           %% test di arresto, tolleranza e maxiter
    i=i+1;
    r=b-A*x;                
    x=x+alpha*(r);                                 %% applico il metodo di richardson
    r=norm(r);
    res(i+1)=(r);                                  %% salvo il residuo
end