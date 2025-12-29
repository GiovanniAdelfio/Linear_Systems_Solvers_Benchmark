function [Q,R] = Fatt_QR(A)
%% Prende in input una matrice A quadrata non singolare, e dà in output una 
%% fattorizzazione QR trovata con le matrici di Householder.
n = size(A,1);
Q = eye(n);
R = A;                                           % Inizializzo R=A e Q come Identità
for i = 1:n-1
    x = R(i:n,i);                                % estraggo il vettore su cui costruire la riflessione
    sig = -sign(x(1))*norm(x,2);                 % scelgo il segno in modo da evitare errori di annullamento
    v = x - sig * eye(n-i+1,1);                  % costruisco v in modo che H*x // e1

    beta = 2 / (v' * v);
    Hi = eye(n);
    Hi(i:n, i:n) = Hi(i:n, i:n) - beta*(v * v'); % costruisco Hi da v

    Q = Q*Hi;                                     % aggiorno Q 
    R = Hi*R;                                     % aggiorno R
end
R = triu(R);                                     %  forza la triangolarità superiore esatta