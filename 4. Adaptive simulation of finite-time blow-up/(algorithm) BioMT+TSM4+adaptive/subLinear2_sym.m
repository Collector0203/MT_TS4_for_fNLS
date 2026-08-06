% This function returns the numerical solution of 1i*\partial_t \psi(x,y,t)
% = (-\Delta)^s \psi(x,y,t) computed by BioMT-Galerkin method, where (x,y)
% \in R^2, t \in R, s \in (0,1).
% 
% input: 
% 1. the values of \psi(x,y,0) on MT-nodes "psi".
% 2. the eigenvectors "E1" of sub-MT-derivative matrix S1.
% 3. truncation parameter "N".
% 4. common matrix "BioMTmatrix_sym", "opt_Sign_SignM", "opt_iTheta",
% "opt_expm_w_sym".
% 
% output: the values of \psi(x,y,dt) on MT-nodes "Psi".
%
% (This program is specifically designed for equation2D_blow_gauss.m.)

function Psi = subLinear2_sym(psi,E1,N,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_expm_w_sym)
MTcoeff_0 = MTcoeff2(psi,N,opt_Sign_SignM,opt_iTheta); % the MT-coefficients of \psi(x,y,0) under {MT(m,mu*x)*MT(n,mu*y)}
Psi = fastPsi_sym(E1,N,MTcoeff_0,opt_expm_w_sym,BioMTmatrix_sym); % using the symmetry of matrices to fast compute "Psi"