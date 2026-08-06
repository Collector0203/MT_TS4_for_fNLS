% This program uses fftn to compute the MT-coefficients
% {a(p,q,r)}_{p,q,r=-N:N-1} of 3D function u(x,y,z), where a(p,q,r) :=
% int(x=-inf:inf)int(y=-inf:inf)int(z=-inf:inf)
% u(x,y,z)*conj(MT(p,x))*conj(MT(q,y))*conj(MT(r,z)).
% 
% input: the function value array of u(x,y,z) on MT-nodes "U", truncation
% "N".
% output: MT-coefficients of u(x,y,z).
%
% (remark): When "U" is computed based on the scaled MT-nodes, "coeff =
% MTcoeff3_origin(U,N)" implies that "u(x,y,z) \approx sum(j,k,l=-N:N-1)
% coeff(j,k,l)*MT(j,mu*x)*MT(k,mu*y)*MT(l,mu*z)".

function coeff = MTcoeff3_origin(U,N)
theta = -pi + pi/N*(0:2*N-1); % generate theta sequence
[theta1,theta2,theta3] = ndgrid(theta,theta,theta);
U1 = U.*(1-1i*tan(theta1/2)).*(1-1i*tan(theta2/2)).*(1-1i*tan(theta3/2));
U2 = fftsum3(U1,N);
[M1,M2,M3] = ndgrid(-N:N-1,-N:N-1,-N:N-1);
coeff = (-1i).^(M1+M2+M3)*pi*sqrt(pi)/16/sqrt(2)/N^3.*U2;