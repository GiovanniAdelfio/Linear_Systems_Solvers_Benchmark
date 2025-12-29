function [L,U,P]= Fatt_LU_piv(M, Tol)
% prende in input una matrice A non singolare, e dà in output una fattorizzazione LU
% con pivoting totale, cioè scelgo sempre il pivot massimo.
% Tol è il valore sotto il quale i pivot, in valore assoluto, sono considerati nulli

U=M;                          
n=size(U,1);
L=eye(n);
P=eye(n);
for j = 1:n
     [mass,k]=max(abs(U(j:n,j)));      %% trovo il pivot massimo
     if k >1
         k=j+k-1;
         if mass < Tol                 %% se il pivot massimo è 'nullo' A è singolare
             error('A è singolare')
         end
         p= 1:n;
         p([j,k])= [k,j];
    
         U=U(p,:);                     %% effettuo le permutazioni su U, P, L
         P=P(p,:);
         L(:,1:j-1)=L(p,1:j-1);
     end

     m=U(j+1:n,j)/U(j,j);              %% trovo il vettore di moltiplicatori 
     L(j+1:n,j)=m;                     %% li salvo in L
     U(j+1:n,:)=U(j+1:n,:)-m*U(j,:);   %% aggiorno U
end 
