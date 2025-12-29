function [x, res] = Gradiente(A,b, Tol)
%% Risolutore di Sistemi Lineari per mezzo del metodo del gradiente
x = zeros(size(b));
r = b-A*x;                              %% residuo iniziale, antigradiente
i = 0;
maxiter = 1000;
res = zeros(maxiter+1,1);
res(1) =  norm(r);
while res(i+1)>Tol && i<maxiter         %% test di arresto
    i = i+1;
    Ar = A*r;                           %% salvo A * r che userò più volte
    alpha = (r'*r)/(Ar'*r);             %% alpha ottimale che minimizza la loss
    x = x + alpha * r;                  %% aggiorno x con il met del gradiente
    r = r - alpha * Ar;                 %% aggiorno il residuo
    res(i+1) = norm(r);                 %% salvo la norma del residuo
end
res = res(1:i+1);                       %% tengo solo i valori non nulli