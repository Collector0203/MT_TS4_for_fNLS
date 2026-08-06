% Considering the 1D nonlinear fractional Schrödinger equation:
% 1i*\partial_t \psi(x,t) = 1/2*(-\Delta)^s \psi(x,t) +
% gamma*|\psi(x,t)|^{2*p}*\psi(x,t), where x \in R, t>0. This program
% solves this equation with parameters (s,p,gamma) = (0.6,1,-1) and initial
% value function \psi_0(x) = exp(-x^2) by BioMT + TSM4, and generates the
% temporal errors about the time step dt.

clear

%% Set the parameters for equations, space discretization, time discretization
s = 0.6; p = 1; gamma = -1; % equation parameter
N = 700; % truncation, \psi(x,t) = sum(n=-N:N-1) a(n,t)*MT(n,mu*x)
mu = 1/sqrt(N); % the scaling factor for MT-functions, change {MT(n,x)} to {MT(n,mu*x)}
w = [0.33780179798991440851, 0.67560359597982881702, ...
    -0.08780179798991440851, -0.85120719195965763405]; % the weights for 4-order time-splitting method
T = 1; % final time
Dt = T./(round(2.^(3:0.5:6))); % time step, t_{n+1} = t_n + dt

%% Generate the MT-derivative product matrix
ind = -N:-1;
L1 = ind.^2 + (ind+1).^2 +(2*ind+1).^2;
L2 = 4*1i*ind.^2; L2 = L2(2:N);
L4 = -ind.*(ind-1); L4 = L4(3:N);
L3 = conj(L2); L5 = L4;
S1 = diag(L1) + diag(L2,-1) + diag(L3,1) + diag(L4,-2) + diag(L5,2); % derivative product matrix
[E1,D1] = eig(S1,'vector'); % compute the eigenvalues D1 and eigenvectors E1 of S1
D2 = D1; E2 = conj(flip(E1,1));

%% Generate the MT-nodes for "mu"
theta = -pi+pi/N*(0:2*N-1)'; node = 1/2*tan(theta/2);
node_mu = node/mu; % MT-nodes for "mu" in 1D (colunm vector)

%% Pre-compute the common matrix to improve the efficiency of codes
% opt_Sign_SignM, used in MTcoeff.m
opt_Sign_SignM = (-1i).^(-N:N-1).'/2/N*sqrt(pi/2).*(-1).^((0:2*N-1).'+mod(N,2));

% opt_iTheta, used in MTcoeff.m
opt_iTheta = 1-1i*tan(theta/2);

% BioMTmatrix, used to compute the value of sum(n=-N:N-1) U(n)*BioMT(n,mu*x) on MT-nodes based on U.
[NN,nnode] = meshgrid(-N:N-1,node);
MTmatrix = MT(NN,nnode);
Mtemp1 = MTmatrix(1:N,1:N);
Mtemp2 = MTmatrix(1:N,N+1:2*N);
Mtemp3 = MTmatrix(N+1:2*N,1:N);
Mtemp4 = MTmatrix(N+1:2*N,N+1:2*N);
BioMTmatrix = [Mtemp1*conj(E1), Mtemp2*conj(E2); Mtemp3*conj(E1), Mtemp4*conj(E2)];

%% Generate the test points
M = 500; pt = -M:M/2000:M; % test interval [-M,M] and test points pt
[NN,ppt] = meshgrid(-N:N-1,pt); MTmatrix = MT(NN,mu*ppt); % common matrix, which is used to compute the value of sum(n=-N:N-1) U(n)*MT(n,mu*x) on test points pt based on U

%% Compute the numerical solution of \psi(x,T)
Unum_coeff = cell(length(Dt),1); % Unum_coeff{k,1} denotes the MT-coefficients of \psi(x,T) corresponding to dt = Dt(k)
Unum_valpt = zeros(length(pt),length(Dt)); % Unum_valpt(:,k) denotes the values of the numerical solution \psi(x,T) corresponding to Dt(k) on test points pt

h = waitbar(0,'please wait ...'); flag = 0; 
for dt = Dt
    %% Generate the initial value function
    U0 = initial(node_mu); % the values of initial function on MT-nodes

    %% Generate the common matrix "opt_expm_w" (used in subLinear.m)
    opt_expm_w2 = exp(-1i*mu^(2*s)*[D1;D2].^s*w(2)*dt);
    opt_expm_w4 = exp(-1i*mu^(2*s)*[D1;D2].^s*w(4)*dt);

    %% Time evolution
    Tstep = round(T/dt); % Tstep represents the number of time steps required from t=0 to t=T
    tic
    for m = 1:Tstep
        Un = U0;
        U1 = subNonlinear(Un,w(1)*dt,gamma,p); % solve the nonlinear subproblem
        U2 = subLinear(U1,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2); % solve the linear subproblem
        U3 = subNonlinear(U2,w(3)*dt,gamma,p);
        U4 = subLinear(U3,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
        U5 = subNonlinear(U4,w(3)*dt,gamma,p);
        U6 = subLinear(U5,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
        U0 = subNonlinear(U6,w(1)*dt,gamma,p);
    end
    toc

    %% Compute and store the numerical solution
    flag = flag+1; waitbar(flag/length(Dt));
    unum_coeff = MTcoeff(U0,N,opt_Sign_SignM,opt_iTheta);
    Unum_coeff{flag,1} = unum_coeff; % MT-coefficients
    unum_nodept = MTmatrix*unum_coeff;
    Unum_valpt(:,flag) = unum_nodept; % the values of numerical solution on test points
end
close(h)

%% Load the reference solution (note the time T in this program and TT in 'ref_solt.mat')
load('ref_solt.mat');
k_ref = -1;
for k = 1:length(TT)
    if abs(TT(k)-T) < 1e-12, k_ref = k; end % find TT(k_ref) that equals to T
end
if k_ref == -1
    disp('error: There is no corresponding reference solution for the termination time T in this program.')
else
    unum_coeff_ref_T = unum_coeff_ref(:,k_ref); % unum_coeff_ref(:,k) denotes the MT-coefficients of \psi_{ref}(x,TT(k)).
end

%% Compute the Linf-error
[NN,ppt] = meshgrid(-N_ref:N_ref-1,pt); MTmatrix = MT(NN,mu_ref*ppt);
unum_nodept_ref = MTmatrix*unum_coeff_ref_T;
error_Linf_T = zeros(length(Dt),1);
for k = 1:length(Dt)
    error_Linf_T(k) = max(abs(unum_nodept_ref-Unum_valpt(:,k)),[],'all');
end

%% Compute the L2-error
theta_ref = -pi+pi/N_ref*(0:2*N_ref-1); 
node_ref = 1/2/mu_ref*tan(theta_ref.'/2); % MT-nodes of N_ref and mu_ref
error_L2_T = zeros(length(Dt),1);
[NN,nnode_ref] = meshgrid(-N:N-1,node_ref); 
for k = 1:length(Dt)
    coeff = Unum_coeff{k,1}; % extract the parameters and coefficients
    uN = MT(NN,mu*nnode_ref)*coeff; % the values of uN(x,T) on MT-nodes corresponding to N_ref
    unum_coeff_T = MTcoeff_origin(uN,N_ref); % the MT-coefficients of uN(x,T) under basis {MT(n,mu_ref*x)}_{-N_ref:N_ref-1}
    error_L2_T(k) = sqrt(abs(1/mu_ref*sum(abs(unum_coeff_T-unum_coeff_ref_T).^2)));
end

%% Roughly plot the error-decay figure
figure
loglog(Dt,error_Linf_T,'*'), hold on, 
loglog(Dt,error_L2_T,'*'), hold on,
loglog(Dt,Dt.^4,'--'),
grid on, xlabel dt, ylabel error,
legend('Linf-error','L2-error','O(dt^{4})')
title(['The numerical time-error of \psi(x,y,',num2str(T),') with N =',num2str(N)])

%% Save the temporal error data
save("error_time.mat","error_Linf_T","error_L2_T","s","p","gamma","N","mu","T","Dt");