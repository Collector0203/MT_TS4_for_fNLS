% This program is specifically designed for equation2D_blow_gauss.m. Its
% effect is equivalent to:
% " BioMTcoeff_0 = [(E1.')*MTcoeff_0(1:N,1:N)*E1, ...
%     (E1.')*MTcoeff_0(1:N,N+1:2*N)*E2; ...
%     (E2.')*MTcoeff_0(N+1:2*N,1:N)*E1, ...
%     (E2.')*MTcoeff_0(N+1:2*N,N+1:2*N)*E2];
% BioMTcoeff_dt = opt_expm_w.*BioMTcoeff_0;
% Psi = BioMTmatrix*BioMTcoeff_dt*(BioMTmatrix.'); ", where
% 1. E2 = conj(flip(E1,1)). 
% 2. opt_tenD_s_sym = opt_tenD_s(1:N,1:N).
% 3. BioMTmatrix_sym = 2*imag(BioMTmatrix([1,N+1:2*N],1:N)).

function Psi = fastPsi_sym(E1,N,MTcoeff_0,opt_expm_w_sym,BioMTmatrix_sym)
MTcoeff_0_sym = MTcoeff_0(1:N,1:N);
BioMTcoeff_0_sym = E1.'*MTcoeff_0_sym*E1;
BioMTcoeff_dt_sym = opt_expm_w_sym.*BioMTcoeff_0_sym;
Psi_temp = -BioMTmatrix_sym*BioMTcoeff_dt_sym*(BioMTmatrix_sym.');
Psi_temp_reo = Psi_temp([1,N+1:-1:3,2:N+1],[1,N+1:-1:3,2:N+1]);
Psi = 1/2*(Psi_temp_reo+Psi_temp_reo.');