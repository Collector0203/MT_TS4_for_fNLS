% This function returns the numerical solution of 1i*\partial_t \psi(x,t) =
% (-\Delta)^s \psi(x,t) computed by BioMT-Galerkin method, where x \in R, t
% \in R, s \in (0,1).
% 
% input: 
% 1. the values of \psi(x,0) on MT-nodes "psi".
% 2. the eigenvectors "E" of MT-derivative matrix S, where E = [E1,0;0,E2].
% 3. truncation parameter "N".
% 4. common matrix "BioMTmatrix", "opt_Sign_SignM", "opt_iTheta", which are
% determined by truncation parameter "N".
% 5. common matrix "opt_expm_w", which is determind by "N", "dt", "mu",
% "s".
% 
% output: the values of \psi(x,dt) on MT-nodes "Psi".
%
% (This function has been optimized for efficiency.)

function Psi = subLinear(psi,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w)
MTcoeff_0 = MTcoeff(psi,N,opt_Sign_SignM,opt_iTheta); % the MT-coefficients of \psi(x,0) under {MT(n,mu*x)}
BioMTcoeff_0 = [E1.'*MTcoeff_0(1:N); E2.'*MTcoeff_0(N+1:2*N)]; % the coefficients under {BioMT(n,mu*x)}
BioMTcoeff_dt = opt_expm_w.*BioMTcoeff_0; % the MT-coefficients of \psi(x,T) under {BioMT(n,mu*x)}
Psi = BioMTmatrix*BioMTcoeff_dt; % compute the values of \psi(x,dt) on MT-nodes through BioMTcoeff_dt directly