% This program uses the symmetry of "coeff" and "MT(NN, mu*nnode)" to fast
% compute "MT(NN,mu*nnode)*coeff*(MT(NN,mu*nnode).')".
% 
% input: 
% 1. MT-coefficients "coeff" generated in equation2D_blow_gauss.m.
% 2. truncation parameter "N".
% 3. scaling factor "mu".
% 4. common data "NN_sym", "nnode_sym", "ell".
% 
% ouput: Target = MT(NN,mu*nnode)*coeff*(MT(NN,mu*nnode).')
%
% (This program is specifically designed for equation2D_blow_gauss.m.)

function Target = Mulcoeff2_sym(coeff,N,mu,NN_sym,nnode_sym,ell)
MTmatrix = MT(NN_sym, mu*nnode_sym);
MTmatrix_rep = 2*real(MTmatrix); MTmatrix_rep(:,2:2:end) = 2*imag(MTmatrix(:,2:2:end));
coeff_rep = (ell.*coeff(1:N,1:N)).*ell.'; 
Target_temp = MTmatrix_rep*coeff_rep*MTmatrix_rep.';
Target = Target_temp([1,N+1:-1:3,2:N+1],[1,N+1:-1:3,2:N+1]);