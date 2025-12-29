function bool = isDiagDominant(A)
% Controlla se A è diagonalmente dominante su tutte le righe e se esiste
% almeno una riga in cui è a dominanza strettamente diagonale

d = abs(diag(A));                    % vettore dei |a_ii|
S = sum(abs(A), 2) - d;              % somma delle altre colonne per ogni riga
bool = all(d >= S);                  % vero solo se |a_ii| >= S(i) per ogni i
if  ~bool || sum(d>S)<1
    bool = false;
end