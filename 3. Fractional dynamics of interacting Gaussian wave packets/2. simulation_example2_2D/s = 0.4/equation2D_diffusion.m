% This program uses BioMT+TSM4 method to solve the 2D nonlinear fractional
% Schrödinger equation, i.e., 1i*\partial_t \psi(x,t) = 1/2*(-\Delta)^s
% \psi(x,t) + gamma*|\psi(x,t)|^{2*p}*\psi(x,t), where x \in R^2, t>0, and
% plots the evolution of |\psi|.
% 
% 1. parameters: (s,p,gamma) = (0.4,1,1).
% 2. initial: \psi_0(x,y) = exp(-8(x-3)^2-8(y-2)^2) + exp(-(x+4)^2-(y+4)^2)
% + exp(-2.5(x-3)^2-2.5(y+4)^2).

clear, close all

%% Set the parameters for equations, space discretization, time discretization
s = 0.4; p = 1; gamma = 1; % equation parameter
N = 256; % truncation, \psi(x,y,t) = sum(m=-N:N-1)sum(n=-N:N-1) a(m,n,t)*MT(m,mu*x)*MT(n,mu*y)
mu = 1/sqrt(N); % the scaling factor for MT-functions, change {MT(m,x)*MT(n,y)} to {MT(m,mu*x)*MT(n,mu*y)}
dt = 0.05; % time step, t_{n+1} = t_n + dt
w = [0.33780179798991440851, 0.67560359597982881702, ...
    -0.08780179798991440851, -0.85120719195965763405]; % the weights for 4-order time-splitting method
T = 18; % final time

%% Generate the MT-derivative product matrix
% Compared to the complete S, E, D, there are: 
% 1. S = [S1,0;0,S2], where S2 = conj(flip(flip(S1,1),2);
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
Tstep = round(T/dt); % "Tstep" represents the number of time steps required from t=0 to t=T
inter = 3; % "Tstep" must be divisible by "inter"
unum_node = zeros(2*N,2*N,inter+1); unum_node(:,:,1) = U0; % unum_node(:,:,q) denotes the numerical values of \psi(x,y,(q-1)*(Tstep/inter)*dt) on MT-nodes

tic, h = waitbar(0,'please wait ...'); teg = 1;
for m = 1:Tstep, waitbar(m/Tstep)
    Un = U0;
    U1 = subNonlinear(Un,w(1)*dt,gamma,p); % solve the nonlinear subproblem
    U2 = subLinear2(U1,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2); % solve the linear subproblem
    U3 = subNonlinear(U2,w(3)*dt,gamma,p);
    U4 = subLinear2(U3,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
    U5 = subNonlinear(U4,w(3)*dt,gamma,p);
    U6 = subLinear2(U5,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
    U0 = subNonlinear(U6,w(1)*dt,gamma,p);
    if mod(m,Tstep/inter) == 0
        teg = teg + 1;
        unum_node(:,:,teg) = U0;
    end
end
close(h), toc

%% Save the necessary data for displaying the evolution
save("simulation2D_diffusion_s=0.4.mat","s","p","gamma","N","mu","dt","T","unum_node","Tstep","inter","opt_Sign_SignM","opt_iTheta");