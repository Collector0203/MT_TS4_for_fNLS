% Considering the 2D nonlinear fractional Schrödinger equation:
% 1i*\partial_t \psi(x,y,t) = 1/2*(-\Delta)^s \psi(x,y,t) +
% gamma*|\psi(x,y,t)|^{2*p}*\psi(x,y,t), where (x,y) \in R^2, t>0. The
% purpose of this program is to verify the following property when solving
% this equation with equation parameters (s,p,gamma) = (0,3,1,1) and
% initial data 1 or 2, using BioMT+TSM4:
% 1. the conservation of discrete mass. Numerical results show that for the
% same termination time T, the deviation of the discrete mass
% |M(t)-M(0)|/M(0), 0<=t<=T, is proportional to the number of time steps,
% indicating that the deviation of the discrete mass is caused solely by
% rounding errors.
% 2. the controllability of the deviation of discrete energy
% |E(t)-E(0)|/E(0).

clear, close all

%% Set the parameters for equations, space discretization, time discretization
s = 0.3; p = 1; gamma = 1; % equation parameter
N = 256; % truncation, \psi(x,y,t) = sum(m=-N:N-1)sum(n=-N:N-1) a(m,n,t)*MT(m,mu*x)*MT(n,mu*y)
mu = 1/sqrt(N); % the scaling factor for MT-functions, change {MT(m,x)*MT(n,y)} to {MT(m,mu*x)*MT(n,mu*y)}
dt = 0.05; % time step, t_{n+1} = t_n + dt
w = [0.33780179798991440851, 0.67560359597982881702, ...
    -0.08780179798991440851, -0.85120719195965763405]; % the weights for 4-order time-splitting method
T = 25; % final time

%% Generate the MT-derivative product matrix
ind = -N:-1;
L1 = ind.^2 + (ind+1).^2 +(2*ind+1).^2;
L2 = 4*1i*ind.^2; L2 = L2(2:N);
L4 = -ind.*(ind-1); L4 = L4(3:N);
L3 = conj(L2); L5 = L4;
S1 = diag(L1) + diag(L2,-1) + diag(L3,1) + diag(L4,-2) + diag(L5,2); % derivative product matrix
[E1,D1] = eig(S1,'vector'); % compute the eigenvalues D1 and eigenvectors E1 of S1
D2 = D1; E2 = conj(flip(E1,1));

%% Generate the MT-nodes
theta = -pi+pi/N*(0:2*N-1); node = 1/2*tan(theta/2);
node_mu = node/mu; [nodey_mu,nodex_mu] = meshgrid(node_mu); % MT-nodes for "mu" in 2D

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
[NN,nnode] = meshgrid(-N:N-1,node_mu);
MTmatrix = MT(NN,mu*nnode);
Mtemp1 = MTmatrix(1:N,1:N);
Mtemp2 = MTmatrix(1:N,N+1:2*N); 
Mtemp3 = MTmatrix(N+1:2*N,1:N); 
Mtemp4 = MTmatrix(N+1:2*N,N+1:2*N);
BioMTmatrix = [Mtemp1*conj(E1), Mtemp2*conj(E2); Mtemp3*conj(E1), Mtemp4*conj(E2)];

% opt_tenD_s, used in energy2D.m
opt_tenD_s = tenD.^s;

% opt_CosTheta, used in mass2D.m and energy2D.m
opt_CosTheta = cos(theta1/2).*cos(theta2/2);

%% Compute and store the discrete mass and energy in time evolution with initial data 1
U0 = initial1(nodex_mu,nodey_mu); % the values of initial function on MT-nodes

Tstep = round(T/dt); % Tstep represents the number of time steps required from t=0 to t=T
Mass1 = zeros(Tstep+1,1); Mass1(1) = mass2D(U0,N,mu,opt_CosTheta); % Mass(q) denotes the discrete mass of \psi(x,y,(q-1)*dt)
Energy1 = zeros(Tstep+1,1); Energy1(1) = energy2D(U0,E1,E2,N,mu,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta); % Energy(q) denotes the discrete energy of \psi(x,y,(q-1)*dt)

h = waitbar(0,'please wait ...');
for m = 1:Tstep, waitbar(m/Tstep)
    Un = U0;
    U1 = subNonlinear(Un,w(1)*dt,gamma,p);
    U2 = subLinear2(U1,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
    U3 = subNonlinear(U2,w(3)*dt,gamma,p);
    U4 = subLinear2(U3,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
    U5 = subNonlinear(U4,w(3)*dt,gamma,p);
    U6 = subLinear2(U5,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
    U0 = subNonlinear(U6,w(1)*dt,gamma,p);
    Mass1(m+1) = mass2D(U0,N,mu,opt_CosTheta);
    Energy1(m+1) = energy2D(U0,E1,E2,N,mu,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta);
end
close(h)

%% Compute and store the discrete mass and energy in time evolution with initial data 2
U0 = initial2(nodex_mu,nodey_mu); % the values of initial function on MT-nodes

Tstep = round(T/dt); % Tstep represents the number of time steps required from t=0 to t=T
Mass2 = zeros(Tstep+1,1); Mass2(1) = mass2D(U0,N,mu,opt_CosTheta); % Mass(q) denotes the discrete mass of \psi(x,y,(q-1)*dt)
Energy2 = zeros(Tstep+1,1); Energy2(1) = energy2D(U0,E1,E2,N,mu,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta); % Energy(q) denotes the discrete energy of \psi(x,y,(q-1)*dt)

h = waitbar(0,'please wait ...');
for m = 1:Tstep, waitbar(m/Tstep)
    Un = U0;
    U1 = subNonlinear(Un,w(1)*dt,gamma,p);
    U2 = subLinear2(U1,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
    U3 = subNonlinear(U2,w(3)*dt,gamma,p);
    U4 = subLinear2(U3,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
    U5 = subNonlinear(U4,w(3)*dt,gamma,p);
    U6 = subLinear2(U5,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
    U0 = subNonlinear(U6,w(1)*dt,gamma,p);
    Mass2(m+1) = mass2D(U0,N,mu,opt_CosTheta);
    Energy2(m+1) = energy2D(U0,E1,E2,N,mu,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta);
end
close(h)

%% Roughly plot the deviation of the discrete mass
figure
drift_mass1 = abs((Mass1-Mass1(1))./Mass1(1));
plot((0:Tstep)*dt,drift_mass1,'*-'), hold on,
drift_mass2 = abs((Mass2-Mass2(1))./Mass2(1));
plot((0:Tstep)*dt,drift_mass2,'*-'), hold on,
xlabel t, ylabel('discrete mass deviation')
legend('initial1','initial2')

%% Roughly plot the deviation of the discrete energy
figure
drift_energy1 = abs((Energy1-Energy1(1))./Energy1(1));
plot((0:Tstep)*dt,drift_energy1,'*-'), hold on,
drift_energy2 = abs((Energy2-Energy2(1))./Energy2(1));
plot((0:Tstep)*dt,drift_energy2,'*-'), hold on,
xlabel t, ylabel('discrete energy deviation')
legend('initial1','initial2')

%% Save the necessary data for the deviation of discrete mass and energy
save("conservation_(0.3,1,1).mat","Mass1","Mass2","Energy1","Energy2","s","p","gamma","N","mu","dt","T","Tstep");