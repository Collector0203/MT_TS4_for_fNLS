% This program uses fft to compute v(k) = sum(j=0:2*N-1)
% u(j)*exp[-1i*k*theta(j)], theta(j) = -pi+pi/N*j, k = -N:N-1. (This
% program requires u to be a column vector.)

function v2 = fftsum(u,N)
v = fft(u);
vv = [v(N+1:2*N);v(1:N)];
v2 = (-1).^((0:2*N-1).'+mod(N,2)).*vv;