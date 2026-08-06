% This program uses fft2 to compute the MT-coefficients
% {a(m,n)}_{m,n=-N:N-1} of 2D function u(x,y), where a(m,n) :=
% int(x=-inf:inf)int(y=-inf:inf) u(x,y)*conj(MT(m,x))*conj(MT(n,y)).
% 
% input: 
% 1. the function value matrix of u(x,y) on MT-nodes "U".
% 2. truncation parameter "N".
% 3. common matrix "opt_Sign_SignM" and "opt_iTheta", which are determined
% by truncation parameter "N".
% 
% ouput: MT-coefficients {coeff(m,n)}, where coeff(m,n) = a(m-N-1,n-N-1).
%
% (This function has been optimized for efficiency.)

function coeff = MTcoeff2(U,N,opt_Sign_SignM,opt_iTheta)
U1 = U.*opt_iTheta;
v = fft2(U1);
vv = [v(N+1:2*N,N+1:2*N), v(N+1:2*N,1:N); v(1:N,N+1:2*N), v(1:N,1:N)];
coeff = opt_Sign_SignM.*vv;