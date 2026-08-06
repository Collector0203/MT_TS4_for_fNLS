% Considering the 2D nonlinear fractional Schrödinger equation:
% 1i*\partial_t \psi(x,y,t) = 1/2*(-\Delta)^s \psi(x,y,t) +
% gamma*|\psi(x,y,t)|^{2*p}*\psi(x,y,t), where (x,y) \in R^2, t>0. This
% program solves this equation with parameters (s,p,gamma) = (0.6,1,-1) and
% initial value function \psi_0(x,y) = exp(-x^2-y^2) by BioMT+TSM4, and
% generates the spatial errors about the truncation parameter N.

clear

%% Set the parameters for equations, space discretization, time discretization
s = 0.6; p = 1; gamma = -1; % equation parameter
NNN = [4, 5, 6, 8, 10, 13, 16, 20, 25, 32, 40, 50, 64, 80, 101, 128]; % truncation, \psi(x,y,t) = sum(m=-N:N-1)sum(n=-N:N-1) a(m,n,t)*MT(m,mu*x)*MT(n,mu*y)
Mu = 1.3./sqrt(NNN); % the scaling factor for MT-functions, change {MT(m,x)*MT(n,y)} to {MT(m,mu*x)*MT(n,mu*y)}
dt = 0.001; % time step, t_{n+1} = t_n + dt
w = [0.33780179798991440851, 0.67560359597982881702, ...
    -0.08780179798991440851, -0.85120719195965763405]; % the weights for 4-order time-splitting method
T = 1; % final time

%% Compute the numerical solution of \psi(x,y,T)
Unum_coeff = cell(length(NNN),1); % Unum_coeff{k,1} denotes the MT-coefficients of the numerical solution \psi(x,y,T) corresponding to NNN(k) and Mu(k)
M = 1000; % test interval [-M,M]
pt = -M:M/2000:M; [pt2,pt1] = meshgrid(pt); % test points [pt1,pt2]
Unum_valpt = zeros(length(pt),length(pt),length(NNN)); % Unum_valpt(:,:,k) denotes the values of the numerical solution \psi(x,y,T) corresponding to NNN(k) and Mu(k) on test points [pt1,pt2]

h = waitbar(0,'please wait ...'); flag = 0;
for N = NNN
    %% Set scaling factor
    flag = flag+1; waitbar(flag/length(NNN))
    mu = Mu(flag); % the scaling factor

    %% Generate the MT-derivative product matrix
    ind = -N:-1;
    L1 = ind.^2 + (ind+1).^2 +(2*ind+1).^2;
    L2 = 4*1i*ind.^2; L2 = L2(2:N);
    L4 = -ind.*(ind-1); L4 = L4(3:N);
    L3 = conj(L2); L5 = L4;
    S1 = diag(L1) + diag(L2,-1) + diag(L3,1) + diag(L4,-2) + diag(L5,2); % derivative product matrix
    [E1,D1] = eig(S1,'vector'); % compute the eigenvalues D1 and eigenvectors E1 of S1
    D2 = D1; E2 = conj(flip(E1,1));

    %% Generate the MT-nodes and the initial value function
    theta = -pi+pi/N*(0:2*N-1); node = 1/2*tan(theta/2);
    node_mu = node/mu; [nodey_mu,nodex_mu] = meshgrid(node_mu); % MT-nodes for "mu" in 2D
    U0 = initial(nodex_mu,nodey_mu); % the values of initial function on MT-nodes

    %% Pre-compute the common matrix to improve the efficiency of codes
    % opt_Sign_SignM, used in MTcoeff2.m
    sign = (-1).^(0:2*N-1);
    [sign1,sign2] = meshgrid(sign,sign);
    opt_Sign = abs(sign1+sign2)-ones(2*N);

    [M2,M1] = meshgrid(-N:N-1,-N:N-1);
    opt_SignM = (-1i).^(M1+M2)*pi/8/N^2;

    opt_Sign_SignM = opt_SignM.*opt_Sign;

    % opt_iTheta, used in MTcoeff2.m
    [theta2,theta1] = meshgrid(theta,theta);
    opt_iTheta = (1-1i*tan(theta1/2)).*(1-1i*tan(theta2/2));

    % opt_expm_w, used in subLinear2.m
    [Dy,Dx] = meshgrid([D1;D2]); tenD = Dx+Dy;
    opt_expm_w2 = exp(-1i*mu^(2*s)*tenD.^s*w(2)*dt);
    opt_expm_w4 = exp(-1i*mu^(2*s)*tenD.^s*w(4)*dt);

    % BioMTmatrix, used to compute the value of sum(m=-N:N-1)sum(n=-N:N-1) U(m,n)*BioMT(m,mu*x)*BioMT(n,mu*y) on MT-nodes based on U.
    [NN,nnode] = meshgrid(-N:N-1,node);
    MTmatrix = MT(NN,nnode);
    Mtemp1 = MTmatrix(1:N,1:N);
    Mtemp2 = MTmatrix(1:N,N+1:2*N);
    Mtemp3 = MTmatrix(N+1:2*N,1:N);
    Mtemp4 = MTmatrix(N+1:2*N,N+1:2*N);
    BioMTmatrix = [Mtemp1*conj(E1), Mtemp2*conj(E2); Mtemp3*conj(E1), Mtemp4*conj(E2)];

    %% Time evolution
    Tstep = round(T/dt); % Tstep represents the number of time steps required from t=0 to t=T
    tic
    for m = 1:Tstep
        Un = U0;
        U1 = subNonlinear(Un,w(1)*dt,gamma,p); % solve the nonlinear subproblem
        U2 = subLinear2(U1,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2); % solve the linear subproblem
        U3 = subNonlinear(U2,w(3)*dt,gamma,p);
        U4 = subLinear2(U3,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
        U5 = subNonlinear(U4,w(3)*dt,gamma,p);
        U6 = subLinear2(U5,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
        U0 = subNonlinear(U6,w(1)*dt,gamma,p);
    end
    toc

    %% Compute and store the numerical solution
    unum_coeff = MTcoeff2(U0,N,opt_Sign_SignM,opt_iTheta);
    Unum_coeff{flag,1} = unum_coeff; % MT-coefficients
    [NN,ppt] = meshgrid(-N:N-1,pt); MTmatrix = MT(NN,mu*ppt); 
    unum_nodept = MTmatrix*unum_coeff*(MTmatrix.');
    Unum_valpt(:,:,flag) = unum_nodept; % the values of numerical solution on test points
end
close(h)

%% Load the reference solution (note the time T in this program and TT in 'ref_solt.mat')
load('ref_solt.mat');
k_ref = -1;
for k = 1:length(TT)
    if abs(TT(k)-T) < 1e-12, k_ref = k; end % find TT(k_ref) that equals to T
end
if k_ref == -1
    error('There is no corresponding reference solution for the termination time T in this program.')
else
    unum_coeff_ref_T = unum_coeff_ref(:,:,k_ref); % unum_coeff_ref(:,:,k) denotes the MT-coefficients of \psi_{ref}(x,y,TT(k)).
end

%% Compute the Linf-error
[NN,ppt] = meshgrid(-N_ref:N_ref-1,pt); MTmatrix = MT(NN,mu_ref*ppt);
unum_nodept_ref = MTmatrix*unum_coeff_ref_T*(MTmatrix.');
error_Linf_T = zeros(length(NNN),1);
for k = 1:length(NNN)
    error_Linf_T(k) = max(abs(unum_nodept_ref-Unum_valpt(:,:,k)),[],'all');
end

%% Compute the L2-error
theta_ref = -pi+pi/N_ref*(0:2*N_ref-1); 
node_ref = 1/2/mu_ref*tan(theta_ref.'/2); % MT-nodes of N_ref and mu_ref
error_L2_T = zeros(length(NNN),1);
for k = 1:length(NNN)
    mu = Mu(k); N = NNN(k); coeff = Unum_coeff{k,1}; % extract the parameters and coefficients
    [NN,nnode_ref] = meshgrid(-N:N-1,node_ref); MTmatrix = MT(NN,mu*nnode_ref);
    uN = MTmatrix*coeff*(MTmatrix.'); % the values of uN(x,y,T) on MT-nodes corresponding to N_ref
    unum_coeff_T = MTcoeff2_origin(uN,N_ref); % the MT-coefficients of uN(x,y,T) under basis {MT(m,mu_ref*x)*MT(n,mu_ref*y)}_{-N_ref:N_ref-1}
    error_L2_T(k) = sqrt(abs(1/mu_ref^2*sum(abs(unum_coeff_T-unum_coeff_ref_T).^2,'all')));
end

%% Roughly plot the error-decay figure
figure
loglog(NNN,error_Linf_T,'*'), hold on,
loglog(NNN,error_L2_T,'*'), hold on,
loglog(NNN,NNN.^(-3-3*s),'--'), hold on
loglog(NNN,NNN.^(-1.5-3*s),'--'), hold on
legend('Linf-error','L2-error','O(N^{-3-3s})','O(N^{-1.5-3s})')
grid on, xlabel N, ylabel error,
title(['The numerical spatial-error of \psi(x,y,',num2str(T),') with dt =',num2str(dt)])

%% Save the spatial error data
save("error_space.mat","error_Linf_T","error_L2_T","s","p","gamma","NNN","Mu","dt","T");