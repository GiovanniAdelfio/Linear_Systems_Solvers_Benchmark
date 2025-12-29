function [L,U,P]= Fatt_LU(M, Tol)
% prende in input una matrice A non singolare, e dà in output una fattorizzazione
% LU con pivoting parziale: cambio pivot solo se 'nullo', poi scelgo il
% primo non nullo.
% Tol è il valore sotto il quale i pivot, in valore assoluto, sono considerati nulli

U=M;                          
n=size(U,1);
L=eye(n);
P=eye(n);
for j = 1:n
    k=j;
    while abs(U(j,j)) <= Tol
         k= k+1;
         if k==n+1                        %% se non ho trovato un pivot non nullo A singolare
            error('A è singolare')
         end
         if abs(U(k,j)) > Tol
             p= 1:n;
             p([j,k])= [k,j];             %% creo il vettore di permutazione

             U=U(p,:);
             P=P(p,:);
             L(:,1:j-1)=L(p,1:j-1);       %% effettuo le permutazioni su U, P, L
             break
         end
     end
     m=U(j+1:n,j)/U(j,j);                 %% trovo il vettore di moltiplicatori 
     L(j+1:n,j)=m;                        %% li salvo in L
     U(j+1:n,:)=U(j+1:n,:)-m*U(j,:);      %% aggiorno U
end 
