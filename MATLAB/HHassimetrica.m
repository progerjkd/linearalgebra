% assimétrica
clc; clear;
%A = [ 4 6 7 2; 1 2 0 1; -2 0 3 -2; 2 1 -2 -1]
A = rand(4);

disp('Transformação de Householder: (obtemos uma matrix de Hessenberg)')
HH = householder(A)
for i=1:50,
[Q,R] = qr(HH); HH = R*Q;
end

disp('Matriz com dentes (ou não) abaixo da diagonal:')
HH

epsilon = 0.001;
n = size(HH,1);
i = 1;

while (i <= n)
    if ((i < n) && (abs(HH(i+1,i)) > epsilon) )
        % seleciona a matrix 2x2 formada com os dentes e calcula os seus
        % autovalores
        M2 = HH(i:i+1,i:i+1);
        [a1, a2] = eqSegGrau(1, -(M2(1,1)+M2(2,2)), det(M2));
        fprintf('Achei o autovalor %1.4f + %1.4fi\n', real(a1), imag(a1));
        fprintf('Achei o autovalor %1.4f + %1.4fi\n', real(a2), imag(a2));
        i = i + 1;
    else
        fprintf('Achei o autovalor %1.4f\n',HH(i,i));
    end
    i = i + 1;
end

% comparando os valores obtidos do nosso algoritmo com os valores da função
% eig() do MATLAB


eig(A)