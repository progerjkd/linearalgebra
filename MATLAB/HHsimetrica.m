% sim�trica
clc; clear;
%A = [ 4 2 2 1; 2 -3 1 1; 2 1 3 1; 1 1 1 2]
% gera uma matriz sim�trica aleat�ria:
N = 5
A = randi(10,N,N) - 1;
A = A - tril(A,-1) + triu(A,1)'

disp('Transforma��o de Householder: (obtemos uma matrix tridiagonal)')
HH = householder(A)
for i=1:50,
[Q,R] = qr(HH); HH = R*Q;
end
HH
Q
disp('Matriz diagonal de autovalores:')
R


% comparando os valores obtidos do nosso algoritmo com os valores da fun��o
% eig() do MATLAB

eig(A)