% This program uses fft2 to compute the MT coefficients {a(m,n)} of 2D
% function u(x,y), which means u(x,y) \appro sum(m=-N:N-1)sum(n=-N:N-1)
% a(m,n)*MT(m,mu*x)*MT(n,mu*y).
% 
% input: the function value matrix of u(x,y) on MT-nodes "U", truncation
% "N", opt-matrix "opt_Sign_SignM" and "opt_iTheta".
% ouput: MT-coefficients {coeff(m,n)}, where coeff(m,n) = a(m-N-1,n-N-1).
% 
% This function has been optimized for efficiency. See
% MT_TSM_equation2D_bo.m for details on "opt_Sign_SignM" and "opt_iTheta".

function coeff = MTcoeff2(U,N,opt_Sign_SignM,opt_iTheta)
U1 = U.*opt_iTheta;
V = fft2(U1);
VV = [V(N+1:2*N,N+1:2*N), V(N+1:2*N,1:N); V(1:N,N+1:2*N), V(1:N,1:N)];
coeff = opt_Sign_SignM.*VV;