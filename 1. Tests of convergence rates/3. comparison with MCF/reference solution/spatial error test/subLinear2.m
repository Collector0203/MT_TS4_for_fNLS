% This function returns the numerical solution of 1i*\partial_t \psi(x,y,t)
% = (-\Delta)^s \psi(x,y,t) computed by BioMT-Galerkin method, where (x,y)
% \in R^2, t \in R, s \in (0,1).
% 
% input: 
% 1. the values of \psi(x,y,0) on MT-nodes "psi".
% 2. the eigenvectors "E" of MT-derivative matrix S, where E = [E1,0;0,E2].
% 3. truncation parameter "N".
% 4. common matrix "BioMTmatrix", "opt_Sign_SignM", "opt_iTheta", which are
% determined by truncation parameter "N".
% 5. common matrix "opt_expm_w", which is determind by "N", "dt", "mu",
% "s".
% 
% output: the values of \psi(x,y,dt) on MT-nodes "Psi".
%
% (This function has been optimized for efficiency.)

function Psi = subLinear2(psi,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w)
MTcoeff_0 = MTcoeff2(psi,N,opt_Sign_SignM,opt_iTheta); % the MT-coefficients of \psi(x,y,0) under {MT(m,mu*x)*MT(n,mu*y)}
BioMTcoeff_0 = [(E1.')*MTcoeff_0(1:N,1:N)*E1, ...
    (E1.')*MTcoeff_0(1:N,N+1:2*N)*E2; ...
    (E2.')*MTcoeff_0(N+1:2*N,1:N)*E1, ...
    (E2.')*MTcoeff_0(N+1:2*N,N+1:2*N)*E2]; % the coefficients under {BioMT(m,mu*x)*BioMT(n,mu*y)}
BioMTcoeff_dt = opt_expm_w.*BioMTcoeff_0; % the MT-coefficients of \psi(x,y,dt) under {BioMT(m,mu*x)*BioMT(n,mu*y)}
Psi = BioMTmatrix*BioMTcoeff_dt*(BioMTmatrix.'); % compute the values of \psi(x,y,dt) on MT-nodes through BioMTcoeff_dt directly