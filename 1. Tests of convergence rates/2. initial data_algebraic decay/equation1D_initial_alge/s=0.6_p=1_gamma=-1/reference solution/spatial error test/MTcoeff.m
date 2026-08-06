% This program uses fft to compute the MT-coefficients {a(k)}_{k=-N:N-1} of
% 1D function u(x), where a(k) := int(x=-inf:inf) u(x)*conj(MT(k,x)).
% 
% input: 
% 1. the function value sequence of u(x) on MT-nodes "U".
% 2. truncation parameter "N".
% 3. common matrix "opt_Sign_SignM" and "opt_iTheta", which are determined
% by truncation parameter "N".
% 
% output: MT-coefficients {coeff(k)} (colunm vector), where coeff(k) =
% a(k-N-1).
%
% (This function has been optimized for efficiency.)

function coeff = MTcoeff(U,N,opt_Sign_SignM,opt_iTheta)
U1 = U.*opt_iTheta;
v = fft(U1);
vv = [v(N+1:2*N);v(1:N)];
coeff = opt_Sign_SignM.*vv;