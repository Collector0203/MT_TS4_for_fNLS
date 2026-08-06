% This program uses fft to compute the MT-coefficients {a(k)}_{k=-N:N-1} of
% 1D function u(x), where a(k) := int(x=-inf:inf) u(x)*conj(MT(k,x)).
% 
% input: the function value sequence of u(x) on MT-nodes "U", truncation
% "N".
% output: MT-coefficients of u(x) (colunm vector).

function coeff = MTcoeff_origin(U,N)
theta = -pi + pi/N*(0:2*N-1).'; % generte theta sequence
U1 = U.*(1-1i*tan(theta/2));
U2 = fftsum(U1,N);
coeff = (-1i).^(-N:N-1).'/2/N*sqrt(pi/2).*U2;