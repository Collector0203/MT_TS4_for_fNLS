% This program uses fftn to compute the MT-coefficients
% {a(p,q,r)}_{p,q,r=-N:N-1} of 3D function u(x,y,z), where a(p,q,r) :=
% int(x=-inf:inf)int(y=-inf:inf)int(z=-inf:inf)
% u(x,y,z)*conj(MT(p,x))*conj(MT(q,y))*conj(MT(r,z)).
% 
% input: 
% 1. the function value array of u(x,y,z) on MT-nodes "U".
% 2. truncation parameter "N".
% 3. common matrix "opt_Sign_SignM" and "opt_iTheta", which are determined
% by truncation parameter N".
% 
% ouput: MT-coefficients {coeff(j,k,l)}, where coeff(j,k,l) =
% a(j-N-1,k-N-1,l-N-1).
% 
% (This function has been optimized for efficiency.)

function coeff = MTcoeff3(U,N,opt_Sign_SignM,opt_iTheta)
U1 = U.*opt_iTheta;
v = fftn(U1);
vv = cat(3,v(:,:,N+1:2*N),v(:,:,1:N));
w1 = vv(1:N,1:N,:);
w2 = vv(1:N,N+1:2*N,:);
w3 = vv(N+1:2*N,1:N,:);
w4 = vv(N+1:2*N,N+1:2*N,:);
w = cat(1,cat(2,w4,w3),cat(2,w2,w1));
coeff = opt_Sign_SignM.*w;