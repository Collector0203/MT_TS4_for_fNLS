% This function returns the numerical solution of 1i*\partial_t \psi(x,y,t)
% = 1/2*(-\Delta)^s \psi(x,y,t) + \gamma*|\psi(x,y,t)|^{2p}*\psi(x,y,t)
% computed by BioMT+TSM4 method, where (x,y) \in R^2, s \in (0,1).
% 
% input:
% 1. the weights for 4-order time-splitting method.
% 2. the values of \psi(x,y,dt) on MT-nodes "Un", generated in
% equation2D_blow_gauss.m.
% 3. equation parameters (s,p,gamma).
% 4. discrete parameters (N,mu,dt).
% 5. the eigenvector matrix "E1" of sub-MT-derivative matrix S1.
% 6. common matrix "BioMTmatrix_sym", "opt_Sign_SignM", "opt_iTheta",
% "opt_tenD_s_sym".
% 
% output: the values of \psi(x,y,dt) on MT-nodes "Un_plus".
%
% (This program is specifically designed for equation2D_blow_gauss.m.)

function Un_plus = BMTS2_sym(w,Un,s,p,gamma,N,mu,dt,E1,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_tenD_s_sym)
opt_expm_w2_sym = exp(-1i*mu^(2*s)*opt_tenD_s_sym*w(2)*dt);
opt_expm_w4_sym = exp(-1i*mu^(2*s)*opt_tenD_s_sym*w(4)*dt);
U1 = subNonlinear(Un,w(1)*dt,gamma,p); % solve the nonlinear subproblem
U2 = subLinear2_sym(U1,E1,N,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_expm_w2_sym); % solve the linear subproblem
U3 = subNonlinear(U2,w(3)*dt,gamma,p);
U4 = subLinear2_sym(U3,E1,N,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_expm_w4_sym);
U5 = subNonlinear(U4,w(3)*dt,gamma,p);
U6 = subLinear2_sym(U5,E1,N,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_expm_w2_sym);
Un_plus = subNonlinear(U6,w(1)*dt,gamma,p);