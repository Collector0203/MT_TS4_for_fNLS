% This program uses fftn to compute v(p,q,r) =
% sum(j=0:2*N-1)sum(k=0:2*N-1)sum(l=0:2*N-1)
% u(j,k,l)*exp[-1i*p*theta(j)]*exp[-1i*q*theta(k)]*exp[-1i*r*theta(l)],
% where theta(j) = -pi+pi/N*j, p,q,r = -N:N-1.

function v2 = fftsum3(u,N)
v = fftn(u);
vv = cat(3,v(:,:,N+1:2*N),v(:,:,1:N));
w1 = vv(1:N,1:N,:); 
w2 = vv(1:N,N+1:2*N,:); 
w3 = vv(N+1:2*N,1:N,:); 
w4 = vv(N+1:2*N,N+1:2*N,:);
w = cat(1,cat(2,w4,w3),cat(2,w2,w1));
sign = (-1).^(0:2*N-1); [sign1,sign2,sign3] = ndgrid(sign,sign,sign); 
v2 = w.*sign1.*sign2.*sign3.*(-1)^mod(N,2);