function D=qr2(A,epsilon)

%Input - A is a symmetric tridiagonal nxn matrix
%      - epsilon is the tolerance
%Output - D is the nx1 vector of eigenvalues

% NUMERICAL METHODS: MATLAB Programs
%(c) 1999 by John H. Mathews and Kurtis D. Fink
%To accompany the textbook:
%NUMERICAL METHODS Using MATLAB,
%by John H. Mathews and Kurtis D. Fink
%ISBN 0-13-270042-5, (c) 1999
%PRENTICE HALL, INC.
%Upper Saddle River, NJ 07458

%Initialize parameters

[n,n]=size(A);
m=n;
D=zeros(n,1);
B=A;

while (m>1)
    m
   while (abs(B(m,m-1))>=epsilon)
     
      %Calculate shift
      disp('B(m-1:m,m-1:m)')
      B(m-1:m,m-1:m)
      S=eig(B(m-1:m,m-1:m))
      B(m,m)
      abs(B(m,m))
      B(m,m)*[1 1]'
      disp('B(m,m)*[1 1]-S')
      B(m,m)*[1 1]'-S
      abs(B(m,m)*[1 1]'-S)
      [j,k]=min([abs(B(m,m)*[1 1]'-S)])
      
      %QR factorization of B
      [Q,U]=qr(B-S(k)*eye(m));
      
      %Calculate next B
      B=U*Q+S(k)*eye(m)      
   end
   
   %Place mth eigenvalue in A(m,m)
   A(1:m,1:m)=B;
   
   %Repeat process on the m-1 x m-1 submatrix of A
   m=m-1;   
   B=A(1:m,1:m);   
end

D=diag(A);

