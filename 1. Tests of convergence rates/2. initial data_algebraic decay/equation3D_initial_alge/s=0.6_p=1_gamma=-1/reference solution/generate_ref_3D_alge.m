% This program generates the reference solution for the 3D fractional
% nonlinear Schrödinger equation with parameters (s,p,gamma) = (0.6,1,-1)
% and initial value function \psi_0(x,y,z) = 1/(1+x^2+y^2+z^2)^(1+sqrt(2)).
% 1. The algorithm parameters "N_ref" is sufficiently large, "dt_ref" is
% sufficiently small.
% 2. This program stores the MT-coefficients of reference solution
% \psi(x,y,z,T) with TT = [0.5,1].

clear

%% Set the parameters for equations, space discretization, time discretization
s = 0.6; p = 1; gamma = -1; % equation parameter
N = 128; % truncation, \psi(x,y,z,t) = sum(j=-N:N-1)sum(k=-N:N-1)sum(l=-N:N-1) a(j,k,l)*MT(j,mu*x)*MT(k,mu*y)*MT(l,mu*z)
mu = 1.5/sqrt(N); % the scaling factor for MT-functions, change {MT(j,x)*MT(k,y)*MT(l,z)} to {MT(j,mu*x)*MT(k,mu*y)*MT(l,mu*z)}
dt = 0.0005; % time step, t_{n+1} = t_n + dt
w = [0.33780179798991440851, 0.67560359597982881702, ...
    -0.08780179798991440851, -0.85120719195965763405]; % the weights for 4-order time-splitting method
T = 1; % final time
TT = [0.5,1]; % the time that record the reference solution

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
E2 = conj(flip(E1,1)); E = [E1,zeros(N);zeros(N),E2]; D = [D1;D1]; % generate the complete E and D

%% Generate the MT-nodes and the initial value function
theta = -pi+pi/N*(0:2*N-1); node = 1/2*tan(theta/2);
node_mu = node/mu; [nodex_mu,nodey_mu,nodez_mu] = ndgrid(node_mu,node_mu,node_mu); % MT-nodes for "mu" in 3D
U0 = initial(nodex_mu,nodey_mu,nodez_mu); % the values of initial function on MT-nodes

%% Pre-compute the common matrix to improve the efficiency of codes
% opt_Sign_SignM, used in MTcoeff3.m
sign = (-1).^(0:2*N-1); 
[sign1,sign2,sign3] = ndgrid(sign,sign,sign); 
opt_Sign = sign1.*sign2.*sign3.*(-1)^mod(N,2);

[M1,M2,M3] = ndgrid(-N:N-1,-N:N-1,-N:N-1);
opt_SignM = (-1i).^(M1+M2+M3)*pi*sqrt(pi)/16/sqrt(2)/N^3;

opt_Sign_SignM = opt_SignM.*opt_Sign; 

% opt_iTheta, used in MTcoeff3.m
[theta1,theta2,theta3] = ndgrid(theta,theta,theta);
opt_iTheta = (1-1i*tan(theta1/2)).*(1-1i*tan(theta2/2)).*(1-1i*tan(theta3/2));

% opt_expm_w, used in subLinear3.m
[Dx,Dy,Dz] = ndgrid(D,D,D); tenD = Dx+Dy+Dz;
opt_expm_w2 = exp(-1i*mu^(2*s)*tenD.^s*w(2)*dt);
opt_expm_w4 = exp(-1i*mu^(2*s)*tenD.^s*w(4)*dt);

% BioMTmatrix, used to compute the value of sum(j=-N:N-1)sum(k=-N:N-1)sum(l=-N:N-1) U(j,k,l)*BioMT(j,mu*x)*BioMT(k,mu*y)*BioMT(l,mu*z) on MT-nodes based on U.
[NN,nnode] = meshgrid(-N:N-1,node); 
MTmatrix = MT(NN,nnode);
BioMTmatrix = MTmatrix*conj(E);

%% Time evolution
Tstep = round(T/dt); % Tstep represents the number of time steps required from t=0 to t=T
unum_node = zeros(2*N,2*N,2*N,length(TT)); % unum_node(:,:,:,q) denotes the numerical values of \psi(x,y,z,(q-1)*(Tstep/inter)*dt) on MT-nodes

tic, h = waitbar(0,'please wait ...');
for m = 1:Tstep, waitbar(m/Tstep)
    Un = U0;
    U1 = subNonlinear(Un,w(1)*dt,gamma,p); % solve the nonlinear subproblem
    U2 = subLinear3(U1,E,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2); % solve the linear subproblem
    U3 = subNonlinear(U2,w(3)*dt,gamma,p);
    U4 = subLinear3(U3,E,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
    U5 = subNonlinear(U4,w(3)*dt,gamma,p);
    U6 = subLinear3(U5,E,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
    U0 = subNonlinear(U6,w(1)*dt,gamma,p);
    for k = 1:length(TT)
        if abs(m*dt-TT(k))<1e-12
            unum_node(:,:,:,k) = U0; % store the reference solution when t = TT(k)
        end
    end
end
close(h), toc

%% Compute the values of numerical solution ​​at test points
mu_ref = mu;
N_ref = N;
dt_ref = dt;
unum_coeff_ref = zeros(2*N,2*N,2*N,length(TT)); % unum_coeff_ref(:,:,:,q) denotes the MT-coefficients of \psi(x,y,z,t=TT(q))
for k = 1:length(TT)
    unum_coeff_ref(:,:,:,k) = MTcoeff3(unum_node(:,:,:,k),N,opt_Sign_SignM,opt_iTheta);
end

%% Save the necessary data for the reference solution
save("ref_solt.mat","s","p","gamma","mu_ref","dt_ref","N_ref","TT","unum_coeff_ref");