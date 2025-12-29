function [x] = Triang_Sup(A,b)
%% Risolutore di sistemi triangolari superiori con matrici quadrate
n = size (A,2);
x = zeros(n,1);
if sum(abs(diag(A)) < 1e-16) > 0                     %% controllo che la diagonale non abbia valori 'nulli'
    error('A è singolare o non triang sup')
end
for i = n:-1:1
    x(i) = (b(i) - A(i,i+1:n) * x(i+1:n))./A(i,i);   %% trovo xi a partire dagli xi già calcolati
end
