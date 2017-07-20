% sim�trica
clc; clear;
%A = [ 4 2 2 1; 2 -3 1 1; 2 1 3 1; 1 1 1 2]
% gera uma matriz sim�trica aleat�ria:
%N = 5
%A = randi(10,N,N) - 1;
%A = A - tril(A,-1) + triu(A,1)'

A = [     7.5431    -6.9605     5.3938     6.1428     9.0157 
   -6.9605    -1.3556     8.3287    -4.5111    -5.7096 
    5.3938     8.3287     5.1765    -6.9538     2.5980 
    6.1428    -4.5111    -6.9538     5.3010     9.9611 
    9.0157    -5.7096     2.5980     9.9611     1.7685 
    ]


disp('Transforma��o de Householder: (obtemos uma matrix tridiagonal)')
HH = householder(A)
for i=1:50,
[Q,R] = clgs(HH); HH = R*Q;
end
HH
Q
disp('Matriz diagonal de autovalores:')
R


% comparando os valores obtidos do nosso algoritmo com os valores da fun��o
% eig() do MATLAB

eig(A)