function [x] = Triang_Inf(A,b)
%% Risolutore di sistemi triangolari inferiori con matrici quadrate
n = size (A,2);
x = zeros(n,1);
if sum(abs(diag(A)) < 1e-17) > 0                    %% controllo che la diagonale non abbia valori 'nulli'
    error('A è singolare o non triang sup')
end
for i = 1:n
    x(i) = (b(i) - A(i,1:i-1) * x(1:i-1))./A(i,i);  %% trovo xi a partire dagli xi già calcolati
end