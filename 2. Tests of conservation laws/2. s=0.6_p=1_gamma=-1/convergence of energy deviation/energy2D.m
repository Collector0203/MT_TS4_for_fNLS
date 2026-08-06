% This function computes the discrete energy of u(x,y) by the values of u
% on the MT-nodes.
% 
% input: 
% 1. the function value matrix of u(x,y) on MT-nodes "U".
% 2. the eigenvectors "E" of MT-derivative matrix S, where E = [E1,0;0,E2].
% 3. truncation parameter "N".
% 4. scaling factor "mu".
% 5. equation parameters: "s", "p", "gamma".
% 6. common matrix "opt_Sign_SignM", "opt_iTheta", "opt_CosTheta", which
% are determined by truncation parameter "N".
% 7. common matrix "opt_tenD_s", which is determined by truncation
% parameter "N" and fractional parameter "s".
% 
% output: discrete energy of u(x,y).
% 
% (This function has been optimized for efficiency.)

function energy = energy2D(U,E1,E2,N,mu,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta)

% compute the discrete kinetic energy
MTcoeff_0 = MTcoeff2(U,N,opt_Sign_SignM,opt_iTheta); % the MT-coefficients of \psi(x,y,0) under {MT(m,mu*x)*MT(n,mu*y)}
BioMTcoeff_0 = [(E1.')*MTcoeff_0(1:N,1:N)*E1, ...
    (E1.')*MTcoeff_0(1:N,N+1:2*N)*E2; ...
    (E2.')*MTcoeff_0(N+1:2*N,1:N)*E1, ...
    (E2.')*MTcoeff_0(N+1:2*N,N+1:2*N)*E2]; % the coefficients under {BioMT(m,mu*x)*BioMT(n,mu*y)}
energy_kinetic_temp = sum(opt_tenD_s.*abs(BioMTcoeff_0).^2,'all');
energy_kinetic = 1/2*mu^(2*s-2)*energy_kinetic_temp;

% compute the discrete potential energy
Utemp = (U.^(p+1))./opt_CosTheta;
energy_potential_temp = (pi/4/mu/N)^2*sum(abs(Utemp).^2,'all');
energy_potential = energy_potential_temp*gamma/(p+1);

% compute the total discrete energy
energy = energy_kinetic + energy_potential;