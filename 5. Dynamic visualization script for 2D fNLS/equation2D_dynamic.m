% This program shows the dynamic evolution of 2D nonlinear fractional
% Schrödinger equation, i.e., 1i*\partial_t \psi(x,t) = 1/2*(-\Delta)^s
% \psi(x,t) + gamma*|\psi(x,t)|^{2*p}*\psi(x,t), where x \in R^2, t>0.
% 
% Q1: What numerical method does this program use?
% A1: This program uses BioMT+TSM4 method to solve the target equation.
% 
% Q2: How to use this program?
% A2: Users need to set the following parameters:
% 1. equation parameters "s","p","gamma", and final time "T", where s \in
% (0,1), p>0, gamma = 1 (defocusing fnls) or -1 (focusing fnls), and T>0.
% 2. initial data, which is set in the script "initial.m".
% 3. discretization parameters "N","mu","dt", where N \in \mathbb{N}, mu>0,
% dt>0. NOTE that "T/dt" must be divisible by 5.
% 4. plotting parameters "M" and "P", where M>0 and P \in \mathbb{N}.
%
% Q3: How does this program display numerical solutions?
% A3: During the time integration, this program outputs the numerical
% solution over the domain [-M,M]^2 in real time. Upon reaching the final
% time T, it pauses for 3 seconds, then plots the numerical solutions at
% the times t = (0:1/5:1)*T (this is why "T/dt" must be divisible by 5).
%
% Q4: Are there any other points to clarify or note?
% A4: The following are some supplementary notes:
% 1. This program uses fixed discretization parameters rather than an
% adaptive strategy, and therefore is not suitable for simulating highly
% resolved blow‑up solutions. However, users can still detect blow‑up by
% visible distortion of the numerical solution and obtain a rough estimate
% of the blow‑up time.
% 2. For home laptops, the recommended value for N is 150-300, and the
% recommended value for mu is N^{-1/2}.

clear, close all

%% Set the parameters for equation, discretization, and plotting
s = 0.7; p = 1; gamma = -1; % equation parameters
T = 10; % final time, solve the equation on (0,T]
N = 150; % discretization-truncation, \psi(x,y,t) = sum(m=-N:N-1)sum(n=-N:N-1) a(m,n,t)*MT(m,mu*x)*MT(n,mu*y)
mu = 1/sqrt(N); % discretization-scaling factor, change {MT(m,x)*MT(n,y)} to {MT(m,mu*x)*MT(n,mu*y)}
dt = 0.05; % discretization-time step, t_{n+1} = t_n + dt
M = 10; P = 400; % display the numerical solution at (P+1)^2 equally spaced nodes on [-M,M]^2

%% Generate the weights for 4-order time-splitting method
w = [0.33780179798991440851, 0.67560359597982881702, -0.08780179798991440851, -0.85120719195965763405];

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

%% Generate the display window for dynamic images
% Display points for dynamic numerical results
pt = -M:2*M/P:M; [pt2,pt1] = meshgrid(pt); % points
[NN,ppt] = meshgrid(-N:N-1,pt); MTmatrix = MT(NN,mu*ppt); % common matrix, which is used to compute the value of sum(m=-N:N-1)sum(n=-N:N-1) U(m,n)*MT(m,mu*x)*MT(n,mu*y) on test points [pt1,pt2] based on U.

% Generate the display window for dynamic images
figure('WindowState', 'maximized');
U0_coeff = MTcoeff2(U0,N,opt_Sign_SignM,opt_iTheta);
U0_nodept = MTmatrix*U0_coeff*(MTmatrix.');

subplot(1,2,1);
h1 = mesh(pt1, pt2, abs(U0_nodept)); % create a mesh object, save its handle as h1
xlabel x, ylabel y, zlabel |\psi(x,y,t)|, grid on,
title('General View','FontSize',18,'Color','k');
colorbar('southoutside'); % place the colorbar below the image

subplot(1,2,2);
h2 = mesh(pt1, pt2, abs(U0_nodept)); % create a mesh object, save its handle as h2
xlabel x, ylabel y,
title('Top View','FontSize',18,'Color','k');
view(2); % set top view
colorbar('southoutside'); % place the colorbar below the image

%% Time evolution (Dynamic visualization)
Tstep = round(T/dt); % Tstep represents the number of time steps required from t=0 to t=T
digit_dt = decimal_num(dt); % the number of decimal places in dt (used in plotting)
unum_node = zeros(2*N,2*N,6); teg = 1; % store the data used to plot the static numerical results
unum_node(:,:,1) = U0; % unum_node(:,:,q) denotes the numerical values of \psi(x,y,(q-1)*(Tstep/5)*dt) on MT-nodes

tic
for m = 1:Tstep
    %% Time evolution from t = (m-1)*Tstep to t = m*Tstep
    Un = U0;
    U1 = subNonlinear(Un,w(1)*dt,gamma,p); % solve the nonlinear subproblem
    U2 = subLinear2(U1,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2); % solve the linear subproblem
    U3 = subNonlinear(U2,w(3)*dt,gamma,p);
    U4 = subLinear2(U3,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w4);
    U5 = subNonlinear(U4,w(3)*dt,gamma,p);
    U6 = subLinear2(U5,E1,E2,N,BioMTmatrix,opt_Sign_SignM,opt_iTheta,opt_expm_w2);
    U0 = subNonlinear(U6,w(1)*dt,gamma,p);

    %% Update the dynamic image
    U0_coeff = MTcoeff2(U0,N,opt_Sign_SignM,opt_iTheta);
    U0_nodept = MTmatrix*U0_coeff*(MTmatrix.');
    U0_nodept_abs = abs(U0_nodept).^2;

    % update the data for subfigure 1
    h1.ZData = U0_nodept_abs; % update the numerical solution 
    h1.CData = U0_nodept_abs; % update color

    % update the data for subfigure 2
    h2.ZData = U0_nodept_abs; % update the numerical solution 
    h2.CData = U0_nodept_abs; % update color

    % update the range of the z-label
    current_max = max(U0_nodept_abs(:));
    margin = current_max*0.05;
    if current_max>1
        new_zlim = [0, current_max+margin];
    else
        new_zlim = [0, 1];
    end
    zlim(subplot(1,2,1), new_zlim);
    zlim(subplot(1,2,2), new_zlim);

    % update the range of the colorbar
    new_clim = [0, current_max+margin];
    clim(subplot(1,2,1), new_clim);
    clim(subplot(1,2,2), new_clim);

    % update the title
    time_str = num2str(m*dt, ['%.', num2str(digit_dt), 'f']);
    sgtitle(['The dynamic evolution of solution \psi(x,y,t), t = ', time_str, ' (final time T = ', num2str(T), ')'], ...
        'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');

    % refresh image
    pause(0.01); % control the refresh rate
    drawnow limitrate;

    %% Stores the data used to plot static numerical results
    if mod(m,Tstep/5) == 0
        teg = teg + 1;
        unum_node(:,:,teg) = U0;
    end
end
toc
pause(3); close;

%% Plot the static numerical results
unum_coeff = zeros(2*N,2*N,6); % unum_coeff(:,:,q) denotes the MT-coefficients of \psi(x,y,(q-1)*(Tstep/5)*dt) under {MT(m,mu*x)*MT(n,mu*y)}
for m = 1:6
    unum_coeff(:,:,m) = MTcoeff2(unum_node(:,:,m),N,opt_Sign_SignM,opt_iTheta);
end
unum_nodept = pagemtimes(pagemtimes(MTmatrix,unum_coeff),MTmatrix.'); % unum_nodept(:,:,q) denotes the values of \psi(x,y,(q-1)*(Tstep/5)*dt) on points [pt1,pt2]

% General view of the numerical solution
figure('WindowState', 'maximized');
Maxunum = max(abs(unum_nodept).^2,[],'all');
for m = 1:6
    subplot(2,3,m)
    mesh(pt1,pt2,abs(unum_nodept(:,:,m)).^2), colorbar
    xlabel x, ylabel y, zlabel |\psi(x,y,t)|, title(sprintf('t = %f', (m-1)*(Tstep/5)*dt));
    axis([-M,M,-M,M,0,1.1*Maxunum]), daspect([1, 1, Maxunum/M/2]);
end
sgtitle('General View', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

% Front view of the numerical solution
figure('WindowState', 'maximized');
Maxunum = max(abs(unum_nodept).^2,[],'all');
for m = 1:6
    subplot(2,3,m)
    mesh(pt1,pt2,abs(unum_nodept(:,:,m)).^2), colorbar
    xlabel x, ylabel y, zlabel |\psi(x,y,t)|, title(sprintf('t = %f', (m-1)*(Tstep/5)*dt));
    axis([-M,M,-M,M,0,1.1*Maxunum]), daspect([1, 1, Maxunum/M/2]), view(0,0);
end
sgtitle('Front View', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

% Top view of the numerical solution
figure('WindowState', 'maximized');
for m = 1:6
    subplot(2,3,m)
    mesh(pt1,pt2,abs(unum_nodept(:,:,m)).^2), colorbar,
    xlabel x, ylabel y, title(sprintf('t = %f', (m-1)*(Tstep/5)*dt));
    axis([-M,M,-M,M]), daspect([1, 1, 1]), view(2);
end
sgtitle('Top View', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');