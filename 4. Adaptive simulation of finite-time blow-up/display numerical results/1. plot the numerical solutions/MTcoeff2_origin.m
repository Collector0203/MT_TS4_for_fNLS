% This program uses fft2 to compute the MT-coefficients
% {a(m,n)}_{m,n=-N:N-1} of 2D function u(x,y), where a(m,n) :=
% int(x=-inf:inf)int(y=-inf:inf) u(x,y)*conj(MT(m,x))*conj(MT(n,y)).
% 
% input: the function value matrix of u(x,y) on MT-nodes "U", truncation
% "N".
% output: MT-coefficients of u(x,y).
%
% (remark): When "U" is computed based on the scaled MT-nodes, "coeff =
% MTcoeff2_origin(U,N)" implies that "u(x,y) \approx sum(j,k=-N:N-1)
% coeff(j,k)*MT(j,mu*x)*MT(k,mu*y)".

function coeff = MTcoeff2_origin(U,N)
theta = -pi + pi/N*(0:2*N-1); % generate theta sequence
[theta2,theta1] = meshgrid(theta,theta);
U1 = U.*(1-1i*tan(theta1/2)).*(1-1i*tan(theta2/2));
U2 = fftsum2(U1,N);
[M2,M1] = meshgrid(-N:N-1,-N:N-1);
coeff = (-1i).^(M1+M2)*pi/8/N^2.*U2;