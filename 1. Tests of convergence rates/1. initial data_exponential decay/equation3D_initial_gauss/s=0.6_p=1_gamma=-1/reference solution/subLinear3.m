% This function returns the numerical solution of 1i*\partial_t
% \psi(x,y,z,t) = (-\Delta)^s \psi(x,y,z,t) computed by BioMT-Galerkin
% method, where (x,y,z) \in R^3, t \in R, s \in (0,1).
% 
% input: 
% 1. the values of \psi(x,y,z,0) on MT-nodes "psi".
% 2. the eigenvectors "E" of MT-derivative matrix S.
% 3. truncation parameter "N".
% 4. common matrix "BioMTmatrix", "opt_Sign_SignM", "opt_iTheta", which are
% determined by truncation parameter "N".
% 5. common matrix "opt_expm_w", which is determind by "N", "dt", "mu",
% "s".
% 
% output: the values of \psi(x,y,z,dt) on MT-nodes "Psi".
% 
% (This function has been optimized for efficiency.)

function Psi = subLinear3(psi,E,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w)
MTcoeff_0 = MTcoeff3(psi,N,opt_Sign_SignM,opt_iTheta); % the MT-coefficients of \psi(x,y,z,0) under {MT(j,mu*x)*MT(k,mu*y)*MT(l,mu*z)}
BioMTcoeff_0_temp = pagemtimes(pagemtimes(E.',MTcoeff_0),E);
BioMTcoeff_0 = tensorprod(BioMTcoeff_0_temp,E.',3,2); % the coefficients under {BioMT(j,mu*x)*BioMT(k,mu*y)*BioMT(l,mu*z)}
BioMTcoeff_dt = opt_expm_w.*BioMTcoeff_0; % the MT-coefficients of \psi(x,y,z,dt) under {BioMT(j,mu*x)*BioMT(k,mu*y)*BioMT(l,mu*z)}
Psi_temp = tensorprod(BioMTcoeff_dt,BioMTmatrix,3,2);
Psi = pagemtimes(pagemtimes(BioMTmatrix,Psi_temp),BioMTmatrix.'); % compute the values of \psi(x,y,z,dt) on MT-nodes through BioMTcoeff_dt directly