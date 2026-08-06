% This program uses fft2 to compute v(n,m) = sum(j=0:2*N-1)sum(k=0:2*N-1)
% u(j,k)*exp[-1i*n*theta(j)]*exp[-1i*m*theta(k)], where theta(j) =
% -pi+pi/N*j, n,m = -N:N-1.

function v2 = fftsum2(u,N)
v = fft2(u);
vv = [v(N+1:2*N,N+1:2*N), v(N+1:2*N,1:N); v(1:N,N+1:2*N), v(1:N,1:N)];
sign = (-1).^(0:2*N-1); [sign1,sign2] = meshgrid(sign,sign); 
Sign = abs(sign1+sign2)-ones(2*N);
v2 = Sign.*vv;