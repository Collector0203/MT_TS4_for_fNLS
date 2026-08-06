% Considering the 2D nonlinear fractional Schrödinger equation:
% 1i*\partial_t \psi(x,y,t) = 1/2*(-\Delta)^s \psi(x,y,t) +
% gamma*|\psi(x,y,t)|^{2*p}*\psi(x,y,t), where (x,y) \in R^2, t>0. The
% purpose of this program is to verify that the deviation of the discrete
% energy satisfies $|E(T)-E(0)|/E(0) = O(dt^4)$ for a fixed final time $T$
% when solving this equation with equation parameters (s,p,gamma) =
% (0.3,1,1) and initial data 1 and 2, using BioMT+TSM4.

clear

%% Set the parameters for equations, space discretization, time discretization
s = 0.3; p = 1; gamma = 1; % equation parameter
N = 256; % truncation, \psi(x,t) = sum(n=-N:N-1) a(n,t)*MT(n,mu*x)
mu = 1/sqrt(N); % the scaling factor for MT-functions, change {MT(n,x)} to {MT(n,mu*x)}
w = [0.33780179798991440851, 0.67560359597982881702, ...
    -0.08780179798991440851, -0.85120719195965763405]; % the weights for 4-order time-splitting method
T = 25; % final time
Dt = 1./round(2.^(4:0.5:6.5)); % time step

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

% BioMTmatrix, used to compute the value of sum(m=-N:N-1)sum(n=-N:N-1) U(m,n)*BioMT(m,mu*x)*BioMT(n,mu*y) on MT-nodes based on U.
[NN,nnode] = meshgrid(-N:N-1,node_mu);
MTmatrix = MT(NN,mu*nnode);
Mtemp1 = MTmatrix(1:N,1:N); 
Mtemp2 = MTmatrix(1:N,N+1:2*N); 
Mtemp3 = MTmatrix(N+1:2*N,1:N); 
Mtemp4 = MTmatrix(N+1:2*N,N+1:2*N);
BioMTmatrix = [Mtemp1*conj(E1), Mtemp2*conj(E2); Mtemp3*conj(E1), Mtemp4*conj(E2)];

% opt_tenD_s, used in energy2D.m
[Dy,Dx] = meshgrid([D1;D2]); tenD = Dx+Dy;
opt_tenD_s = tenD.^s;

% opt_CosTheta, used in energy2D.m
opt_CosTheta = cos(theta1/2).*cos(theta2/2);

%% Vary the time step and record the energy deviation at final time for initial data 1
flag = 0; % time step, t_{n+1} = t_n + dt
energy_drift_1 = zeros(length(Dt),1);

for dt = Dt
    flag = flag + 1;

    %% Generate the initial values ​​and the corresponding descrete energy
    U0 = initial1(nodex_mu,nodey_mu); % the values of initial function on MT-nodes
    Energy_initial = energy2D(U0,E1,E2,N,mu,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta);

    %% Generate the common matrix "opt_expm_w"
    opt_expm_w2 = exp(-1i*mu^(2*s)*tenD.^s*w(2)*dt);
    opt_expm_w4 = exp(-1i*mu^(2*s)*tenD.^s*w(4)*dt);

    %% Generate the descrete energy of the numerical solution at final time
    Tstep = round(T/dt); % Tstep represents the number of time steps required from t=0 to t=T
    h = waitbar(0,['The number of tasks remaining: ',num2str(length(Dt)-flag)]);
    for m = 1:Tstep, waitbar(m/Tstep)
        Un = U0;
        U1 = subNonlinear(Un,w(1)*dt,gamma,p);
        U2 = subLinear2(U1,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
        U3 = subNonlinear(U2,w(3)*dt,gamma,p);
        U4 = subLinear2(U3,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
        U5 = subNonlinear(U4,w(3)*dt,gamma,p);
        U6 = subLinear2(U5,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
        U0 = subNonlinear(U6,w(1)*dt,gamma,p);
    end
    close(h)
    Energy_final = energy2D(U0,E1,E2,N,mu,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta);

    %% Compute the energy deviation
    energy_drift_1(flag) = abs(Energy_final-Energy_initial)/Energy_initial;
end

%% Vary the time step and record the energy deviation at final time for initial data 2
flag = 0; % time step, t_{n+1} = t_n + dt
energy_drift_2 = zeros(length(Dt),1);

for dt = Dt
    flag = flag + 1;

    %% Generate the initial values ​​and the corresponding descrete energy
    U0 = initial2(nodex_mu,nodey_mu); % the values of initial function on MT-nodes
    Energy_initial = energy2D(U0,E1,E2,N,mu,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta);

    %% Generate the common matrix "opt_expm_w"
    opt_expm_w2 = exp(-1i*mu^(2*s)*tenD.^s*w(2)*dt);
    opt_expm_w4 = exp(-1i*mu^(2*s)*tenD.^s*w(4)*dt);

    %% Generate the descrete energy of the numerical solution at final time
    Tstep = round(T/dt); % Tstep represents the number of time steps required from t=0 to t=T
    h = waitbar(0,['The number of tasks remaining: ',num2str(length(Dt)-flag)]);
    for m = 1:Tstep, waitbar(m/Tstep)
        Un = U0;
        U1 = subNonlinear(Un,w(1)*dt,gamma,p);
        U2 = subLinear2(U1,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
        U3 = subNonlinear(U2,w(3)*dt,gamma,p);
        U4 = subLinear2(U3,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
        U5 = subNonlinear(U4,w(3)*dt,gamma,p);
        U6 = subLinear2(U5,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
        U0 = subNonlinear(U6,w(1)*dt,gamma,p);
    end
    close(h)
    Energy_final = energy2D(U0,E1,E2,N,mu,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta);

    %% Compute the energy deviation
    energy_drift_2(flag) = abs(Energy_final-Energy_initial)/Energy_initial;
end

%% Roughly plot the energy deviation at final time T
loglog(Dt,energy_drift_1,'*'), hold on,
loglog(Dt,energy_drift_2,'*'), hold on,
loglog(Dt,10^-2.5*Dt.^4,'--')
xlabel dt, ylabel |E(T)-E(0)|/E(0)

%% Save the necessary data for the deviation of discrete energy
save("energy_order_(0.3,1,1).mat","energy_drift_1","energy_drift_2","s","p","gamma","N","mu","T","Dt");