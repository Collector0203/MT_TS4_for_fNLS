% This function computes the "Findex" of the 2D MT-coefficients, which
% characterizes the proportion of the tail coefficient in "coeff" among all
% coefficients.
%
% input: 2D MT-coefficients "coeff", which size (2*N)*(2*N); truncation
% "N", which means the index -N:N-1.

function Findex = findex2D(coeff,N)
coeff_half = coeff(N+1:2*N,N+1:2*N);
M = round(2/3*N);
f_part1 = sum(abs(coeff_half(M:N,1:M)).^2,'all');
f_part2 = sum(abs(coeff_half(1:M,M:N)).^2,'all');
f_part3 = sum(abs(coeff_half(M+1:N,M+1:N)).^2,'all');
f_total = sum(abs(coeff_half).^2,'all');
Findex = (f_part1 + f_part2 + f_part3)/f_total;