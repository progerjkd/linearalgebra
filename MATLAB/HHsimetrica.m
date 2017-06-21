% sim�trica
clc; clear;
A = [ 4 2 2 1; 2 -3 1 1; 2 1 3 1; 1 1 1 2]
A = [ -2.000000 2.000000 7.000000 -4.000000; 2.000000 8.000000 9.000000 -2.000000; 7.000000 9.000000 -5.000000 5.000000; -4.000000 -2.000000 5.000000 0.000000]
% gera uma matriz sim�trica aleat�ria:
%N = 5
%A = randi(10,N,N) - 1;
%A = A - tril(A,-1) + triu(A,1)'

disp('Transforma��o de Householder: (obtemos uma matrix tridiagonal)')
HH = householder(A)
HH2 = HH;
for i=1:50,
[Q,R] = qr(HH); HH = R*Q;
end
HH
Q
disp('Matriz diagonal de autovalores:')
R



disp('qr2:')
eps = 1e-10
HH2
qr2(HH2, eps)
% comparando os valores obtidos do nosso algoritmo com os valores da fun��o
% eig() do MATLAB


eig(A)