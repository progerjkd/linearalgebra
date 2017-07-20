function [A,d] = qrhous(A);
% QRHOUS   computes the  QR decomposition A = QR
%          using Householder transformations.
%          The output matrix A contains the Householder
%          vectors u and the upper triangle of R. The
%          diagonal of R is stored in vector d.
[m,n]=size(A);
for j = 1:n,
s = norm(A(j:m,j));
if s == 0 , error('rank(A) < n'), end
if A(j,j) >= 0 , d(j)=-s; else d(j)=s; end
fak = sqrt(s*(s+abs(A(j,j))));
A(j,j) = A(j,j)-d(j);
A(j:m,j) = A(j:m,j)/fak;
%   transformation of the rest of
%   the matrix G := G - u*(u’*G)
if j<n,
A(j:m,j+1:n) = A(j:m,j+1:n) - ...
A(j:m,j)*(A(j:m,j)'*A(j:m,j+1:n));
end
end