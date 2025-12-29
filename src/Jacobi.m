function [x,r_jac] = Jacobi(A,b,Tol)
%% Risolutore di sistemi lineari mediante il metodo di Jacobi
i = 0;
maxiter = 1000;
r_jac = zeros(1,maxiter);
x = zeros(size(b));
D = diag(A);                    %% Salvo D come un vettore essendo sparsa
E = A;                          
E(1:size(A,1)+1:end) = zeros(size(A,1),1);  %% creo E=A, e annullo la sua diagonale
r = norm(b-A*x);
r_jac(1) = r; 
while r>Tol && i<maxiter        %% condizioni di arresto, Tolleranza e maxiter
    i=i+1;
    x=(b-E*x)./D;               %% applico il metodo, trovo xi da xi-1
    r=norm(b-A*x);              %% calcolo il nuovo residuo
    r_jac(i+1)=r;               %% salvo il nuovo residuo
end