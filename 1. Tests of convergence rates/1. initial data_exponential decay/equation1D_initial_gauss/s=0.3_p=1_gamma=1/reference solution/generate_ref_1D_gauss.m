% This program generates the reference solution for the 1D fractional
% nonlinear Schrödinger equation with parameters (s,p,gamma) = (0.3,1,1)
% and initial value function \psi_0(x) = exp(-x^2).
% 1. The algorithm parameters "N_ref" is sufficiently large, "dt_ref" is
% sufficiently small.
% 2. This program stores the MT-coefficients of reference solution
% \psi(x,T) with TT = [1,2,3].

clear

%% Set the parameters for equations, space discretization, time discretization
s = 0.3; p = 1; gamma = 1; % equation parameter
N = 700; % truncation, \psi(x,t) = sum(n=-N:N-1) a(n,t)*MT(n,mu*x)
mu = 0.6/sqrt(N); % the scaling factor for MT-functions, change {MT(n,x)} to {MT(n,mu*x)}
dt = 0.0005; % time step, t_{n+1} = t_n + dt
w = [0.33780179798991440851, 0.67560359597982881702, ...
    -0.08780179798991440851, -0.85120719195965763405]; % the weights for 4-order time-splitting method
T = 3; % final time
TT = [1,2,3]; % the time that record the reference solution

%% Generate the MT-derivative product matrix
% Compared to the complete S, E, D, there are:
% 1. S = [S1,0;0,S2], where S2 = conj(flip(flip(S1,1),2));
% 2. E = [E1,0;0,E2], where E2 = conj(flip(E1,1));
% 3. D = [D1;D2], where D2 = D1 and S1*E1 = E1*diag(D1);
ind = -N:-1;
L1 = ind.^2 + (ind+1).^2 +(2*ind+1).^2;
L2 = 4*1i*ind.^2; L2 = L2(2:N);
L4 = -ind.*(ind-1); L4 = L4(3:N);
L3 = conj(L2); L5 = L4;
S1 = diag(L1) + diag(L2,-1) + diag(L3,1) + diag(L4,-2) + diag(L5,2); % derivative product matrix
[E1,D1] = eig(S1,'vector'); % compute the eigenvalues D1 and eigenvectors E1 of S1
D2 = D1; E2 = conj(flip(E1,1));

%% Generate the MT-nodes and the initial value function
theta = -pi+pi/N*(0:2*N-1)'; node = 1/2*tan(theta/2);
node_mu = node/mu; % MT-nodes for "mu" in 1D (colunm vector)
U0 = initial(node_mu); % the values of initial function on MT-nodes

%% Pre-compute the common matrix to improve the efficiency of codes
% opt_Sign_SignM, used in MTcoeff.m
opt_Sign_SignM = (-1i).^(-N:N-1).'/2/N*sqrt(pi/2).*(-1).^((0:2*N-1).'+mod(N,2));

% opt_iTheta, used in MTcoeff.m
opt_iTheta = 1-1i*tan(theta/2);

% opt_expm_w, used in subLinear.m
opt_expm_w2 = exp(-1i*mu^(2*s)*[D1;D2].^s*w(2)*dt);
opt_expm_w4 = exp(-1i*mu^(2*s)*[D1;D2].^s*w(4)*dt);

% BioMTmatrix, used to compute the value of sum(n=-N:N-1) U(n)*BioMT(n,mu*x) on MT-nodes based on U.
[NN,nnode] = meshgrid(-N:N-1,node);
MTmatrix = MT(NN,nnode);
Mtemp1 = MTmatrix(1:N,1:N);
Mtemp2 = MTmatrix(1:N,N+1:2*N);
Mtemp3 = MTmatrix(N+1:2*N,1:N);
Mtemp4 = MTmatrix(N+1:2*N,N+1:2*N);
BioMTmatrix = [Mtemp1*conj(E1), Mtemp2*conj(E2); Mtemp3*conj(E1), Mtemp4*conj(E2)];

%% Time evolution
Tstep = round(T/dt); % Tstep represents the number of time steps required from t=0 to t=T
unum_node = zeros(2*N,length(TT)); % unum_node(:,q) denotes the numerical values of \psi(x,t=TT(q)) on MT-nodes

tic, h = waitbar(0,'please wait ...');
for m = 1:Tstep, waitbar(m/Tstep)
    Un = U0;
    U1 = subNonlinear(Un,w(1)*dt,gamma,p); % solve the nonlinear subproblem
    U2 = subLinear(U1,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2); % solve the linear subproblem
    U3 = subNonlinear(U2,w(3)*dt,gamma,p);
    U4 = subLinear(U3,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
    U5 = subNonlinear(U4,w(3)*dt,gamma,p);
    U6 = subLinear(U5,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
    U0 = subNonlinear(U6,w(1)*dt,gamma,p);
    for k = 1:length(TT)
        if abs(m*dt-TT(k))<1e-12
            unum_node(:,k) = U0; % store the reference solution when t = TT(k)
        end
    end
end
close(h), toc

%% Compute the values of numerical solution ​​at test points
mu_ref = mu;
N_ref = N;
dt_ref = dt;
unum_coeff_ref = zeros(2*N,length(TT)); % unum_coeff_ref(:,q) denotes the MT-coefficients of \psi(x,t=TT(q))
for k = 1:length(TT)
    unum_coeff_ref(:,k) = MTcoeff(unum_node(:,k),N,opt_Sign_SignM,opt_iTheta);
end

%% Save the necessary data for the reference solution
save("ref_solt.mat","s","p","gamma","mu_ref","dt_ref","N_ref","TT","unum_coeff_ref");