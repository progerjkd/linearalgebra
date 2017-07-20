% function [Q,R] = clgs(A)
% % CLGS  classical Gram-Schmidt orthogonalization
% [m,n] = size(A); R = zeros(n);
% Q=A;
% for k = 1:n,
% for i = 1:k-1,
% R(i,k) = Q(:,i)'*Q(:,k);
% end                                 % remove for
% for i = 1:k-1,                      % modified-Gram-Schmidt
% Q(:,k) = Q(:,k)-R(i,k)*Q(:,i);
% end
% R(k,k) = norm(Q(:,k)); Q(:,k) = Q(:,k)/R(k,k);
% end


function [Q,R]=clgs(A)
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n,n);
for j=1:n
V=A(:,j);
for i=1:j-1
R(i,j)=Q(:,i)'*A(:,j);
V=V-R(i,j)*Q(:,i);
end
R(j,j)=norm(V);
Q(:,j)=V/R(j,j);
end