% This program uses BioMT + TSM4 + adaptive method to solve the finite-time
% blow-up problem in 2D fNLS, i.e., 1i*\partial_t \psi(x,t) =
% 1/2*(-\Delta)^s \psi(x,t) + gamma*|\psi(x,t)|^{2*p}*\psi(x,t), where x
% \in R^2, t>0.
% 
% 1. parameters: (s,p,gamma) = (0.7,1.5,-1).
% 2. initial data: \psi_0(x,y) = exp(-x^2/2-y^2/2). 
%
% (This program has been optimized for efficiency by using the real
% symmetry of the initial data.)

clear, close all

%% Set equation parameters and truncation parameter
s = 0.7; p = 1.5; gamma = -1; % equation parameters
N = 256; % truncation parameter
w = [0.33780179798991440851, 0.67560359597982881702, ...
    -0.08780179798991440851, -0.85120719195965763405]; % the weights for 4-order time-splitting method

%% Set adaptive parameters
U_max_bound = 20.5; % (termination condition) time iterates until U_max > U_max_bound
eps_time = 1e-12; % the tolerance error corresponding to the time step during the adaptive process
safe_factor = 0.95; % the safety factor used when adjusting the time step
eps_space = 1e-14; % the tolerance error corresponding to the scaling factor during the adaptive process
ratio = 0.99; % the ratio used to adjust the scaling factor

%% Generate the MT-nodes in 2D
theta = -pi+pi/N*(0:2*N-1); node = 1/2*tan(theta/2);
[nodey,nodex] = meshgrid(node);

%% Generate the MT-derivative product matrix
load('eigS1_N=256.mat');
E2 = conj(flip(E1,1)); D2 = D1; % generate E2 and D2
S = [S1, zeros(N); zeros(N), conj(flip(flip(S1,1),2))]; % generate the complete MT-derivative product matrix

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

% BioMTmatrix_sym, used in BMTS2_sym.m
[NN,nnode] = meshgrid(-N:N-1,node);
MTmatrix = MT(NN,nnode);
Mtemp1 = MTmatrix(1:N,1:N); 
Mtemp2 = MTmatrix(1:N,N+1:2*N); 
Mtemp3 = MTmatrix(N+1:2*N,1:N); 
Mtemp4 = MTmatrix(N+1:2*N,N+1:2*N);
BioMTmatrix = [Mtemp1*conj(E1), Mtemp2*conj(E2); Mtemp3*conj(E1), Mtemp4*conj(E2)];
BioMTmatrix_sym = 2*imag(BioMTmatrix([1,N+1:2*N],1:N));

% opt_tenD_s, used in energy2D.m; opt_tenD_s_sym, used in subLinear2_sym.m, BMTS2_sym.m
[Dy,Dx] = meshgrid([D1;D2]); tenD = Dx+Dy;
opt_tenD_s = tenD.^s;
opt_tenD_s_sym = opt_tenD_s(1:N,1:N);

% opt_CosTheta, used in mass2D.m, energy2D.m
opt_CosTheta = cos(theta1/2).*cos(theta2/2);

% NN_sym, nnode_sym, ell, used in Mulcoeff2_sym.m
idx_row = [1,N+1:2*N]; idx_col = 1:N;
NN_sym = NN(idx_row,idx_col);
nnode_sym = nnode(idx_row,idx_col);
ell = ones(N,1); ell(2:2:end) = 1i;

%% Store discrete and evolution parameters
% time parameters storage ("Tnum" denotes the number of moment)
dt_initial = 0.01; % time step (will adaptive in subsequent)
TT = 0; % TT(q) = t_{q-1}, i.e., TT(1) = 0, TT(2) = Dt(1), TT(3) = Dt(1)+Dt(2), ... (q = 1, ..., Tnum)

% space parameters storage
mu_initial = 1/sqrt(N); % scaling factor (will adaptive in subsequent)
Mu = mu_initial; % Mu(q) denotes the scaling factor corresponding to \psi(x,y,TT(q)) (q = 1, ..., Tnum)

nodex_initial = nodex/mu_initial; nodey_initial = nodey/mu_initial; % MT-nodes of mu
unum_node_initial = initial(nodex_initial,nodey_initial); % the values of initial function on MT-nodes
unum_coeff_initial = MTcoeff2(unum_node_initial,N,opt_Sign_SignM,opt_iTheta); % the MT-coefficients of initial function
Findex = findex2D(unum_coeff_initial,N); % Findex(q) denotes the findex of the MT-coefficients of \psi(x,y,TT(q)) (q = 1, ..., Tnum)

% solution parameters storage
Mass = mass2D(unum_node_initial,N,mu_initial,opt_CosTheta); % Mass(q) denotes the discrete mass of \psi(x,y,TT(q)) (q = 1, ..., Tnum)
Energy = energy2D(unum_node_initial,E1,E2,N,mu_initial,s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta); % Energy(q) denotes the discrete energy of \psi(x,y,TT(q)) (q = 1, ..., Tnum)
U_max = norm(unum_node_initial(:),inf); % U_max(:,q) denotes the Linf-norm of the numerical values of \psi(x,y,TT(q)) (q = 1, ..., Tnum)
Ugrad_L2_temp1 = trace(unum_coeff_initial.'*S*conj(unum_coeff_initial));
Ugrad_L2_temp2 = trace(unum_coeff_initial*S*unum_coeff_initial');
Ugrad_L2 = sqrt(abs(Ugrad_L2_temp1 + Ugrad_L2_temp2)); % Ugrad_L2(:,q) denotes the L2-norm of \nabla \psi(x,y,TT(q)) (q = 1, ..., Tnum)

%% Time evolution
flag = 1; % pedometer

while U_max(flag) < U_max_bound
    %% Read the data of previous step
    mu_pre = Mu(flag); % scaling factor of previous step
    if flag == 1
        unum_node_pre = unum_node_initial; % node values of previous step
        dt_pre = dt_initial;
    else
        dt_pre = TT(flag)-TT(flag-1); % time step of previous step
    end

    %% Adjust the time step based on the scaling factor "mu_pre" from the previous step, in accordance with the "step-doubling"
    t_start = tic; % start timing this iteration

    Un_plus_1 = BMTS2_sym(w,unum_node_pre,s,p,gamma,N,mu_pre,dt_pre,E1,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_tenD_s_sym); % complete step "dt"
    Un_plus_half = BMTS2_sym(w,unum_node_pre,s,p,gamma,N,mu_pre,dt_pre/2,E1,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_tenD_s_sym); % half step "dt"
    Un_plus_2 = BMTS2_sym(w,Un_plus_half,s,p,gamma,N,mu_pre,dt_pre/2,E1,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_tenD_s_sym); % half step "dt"
    error_local = norm((Un_plus_2(:)-Un_plus_1(:)),inf)/15; % error_local = O(dt^5)

    while error_local > eps_time
        dt_now = safe_factor*(eps_time/error_local)^(1/5)*dt_pre; % adjust the time step
        Un_plus_1 = BMTS2_sym(w,unum_node_pre,s,p,gamma,N,mu_pre,dt_now,E1,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_tenD_s_sym); % complete step "dt"
        Un_plus_half = BMTS2_sym(w,unum_node_pre,s,p,gamma,N,mu_pre,dt_now/2,E1,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_tenD_s_sym); % half step "dt"
        Un_plus_2 = BMTS2_sym(w,Un_plus_half,s,p,gamma,N,mu_pre,dt_now/2,E1,BioMTmatrix_sym,opt_Sign_SignM,opt_iTheta,opt_tenD_s_sym); % half step "dt"
        error_local = norm((Un_plus_2(:)-Un_plus_1(:)),inf)/15;
        dt_pre = dt_now;
    end
    dt_now = dt_pre;
    Un_plus_ref = Un_plus_2;

    TT(flag+1) = TT(flag) + dt_now; % store new time parameters

    %% Adjust the scaling factor based on the new time step "dt_now", in accordance with the "tail criterion"
    coeff_ref = MTcoeff2(Un_plus_ref,N,opt_Sign_SignM,opt_iTheta); % the MT-coefficients corresponding to "Un_plus_ref" and "mu_pre"
    findex_ref = findex2D(coeff_ref,N); % findex for reference

    if findex_ref > eps_space % if "findex_ref" exceeds "eps_space", begin the adaptation of scaling factor
        % Three possible outcomes may occur when adjusting the scaling factor; we distinguish them using "teg = 1, 2, 3", respectively:
        % 1. Only a single adjustment step was taken, but the new "findex" was larger than the old one. Therefore, "mu" does not need adjustment.
        % 2. A smaller "findex" was obtained after one adjustment step; adjustment continued until "findex < eps_space".
        % 3. A smaller "findex" was obtained after one adjustment step; adjustment continued until the new "findex" was larger than the "findex" from the previous step.

        mu_inc = mu_pre/ratio; % increase the scaling factor
        unum_node_inc = Mulcoeff2_sym(coeff_ref,N,mu_pre/mu_inc,NN_sym,nnode_sym,ell); % the values of the numerical solution based on "Un_plus_ref" at the MT-nodes for "mu_inc"
        coeff_inc = MTcoeff2(unum_node_inc,N,opt_Sign_SignM,opt_iTheta); 
        findex_inc = findex2D(coeff_inc,N); % new findex
        
        if findex_inc >= findex_ref % first case: no need for adjustment
            Mu(flag+1) = mu_pre; % store new scaling factor
            unum_node_pre = Un_plus_ref; % update the node values for next step
        else
            teg = 0;
            while findex_inc < findex_ref
                findex_ref = findex_inc; % updata "findex_ref"
                mu_inc_pre = mu_inc; unum_node_inc_pre = unum_node_inc; % store "mu" and "node values" from the previous step
                mu_inc_new = mu_inc/ratio; % new "mu"
                unum_node_inc_new = Mulcoeff2_sym(coeff_ref,N,mu_pre/mu_inc_new,NN_sym,nnode_sym,ell); % new "node values"
                mu_inc = mu_inc_new; unum_node_inc = unum_node_inc_new; % updata "mu" and "node values"
                coeff_inc_new = MTcoeff2(unum_node_inc_new,N,opt_Sign_SignM,opt_iTheta); 
                findex_inc = findex2D(coeff_inc_new,N); % new findex
                if findex_inc < eps_space
                    Mu(flag+1) = mu_inc;
                    unum_node_pre = unum_node_inc; % second case: adjustment continued until "findex < eps_space"
                    teg = 1;
                    break;
                end
            end
            if teg == 0 % third case: adjustment continued until the new "findex" was larger than the "findex" from the previous step
                Mu(flag+1) = mu_inc_pre; % go back to the previous scaling factor
                unum_node_pre = unum_node_inc_pre; % the numerical solution at time t_n, with the scaling factor already adjusted
            end
        end
    
    else % if "findex_ref" doesn't exceed "eps_space", disable the adaptation of scaling factor
        Mu(flag+1) = mu_pre; % store new scaling factor
        unum_node_pre = Un_plus_ref; % update the node values for next step
    end

    time_iter(flag) = toc(t_start); % stop timing this iteration, and store the cost time

    %% Store new "Findex", "Mass", "Energy", "U_max", "Ugrad_L2", "time_iter"
    coeff_temp = MTcoeff2(unum_node_pre,N,opt_Sign_SignM,opt_iTheta);
    Findex(flag+1) = findex2D(coeff_temp,N); % store the findex of \psi(x,y,TT(flag+1))
    Mass(flag+1) = mass2D(unum_node_pre,N,Mu(flag+1),opt_CosTheta); % store mass of \psi(x,y,TT(flag+1))
    Energy(flag+1) = energy2D(unum_node_pre,E1,E2,N,Mu(flag+1),s,p,gamma,opt_Sign_SignM,opt_iTheta,opt_tenD_s,opt_CosTheta); % store energy of \psi(x,y,TT(flag+1))
    U_max(flag+1) = norm(unum_node_pre(:),inf); % store the Linf-norm of \psi(x,y,TT(flag+1))
    Ugrad_L2_temp1 = trace(coeff_temp.'*S*conj(coeff_temp));
    Ugrad_L2_temp2 = trace(coeff_temp*S*coeff_temp');
    Ugrad_L2(flag+1) = sqrt(abs(Ugrad_L2_temp1 + Ugrad_L2_temp2)); % store the L2-norm of \nabla \psi(x,y,TT(flag+1))
    flag = flag + 1; % update the pedometer

    %% Store MT-coefficients for plotting
    if flag == 797
        coeff_store{1,1} = flag;
        coeff_store{1,2} = unum_node_pre;
    elseif flag == 1622
        coeff_store{2,1} = flag;
        coeff_store{2,2} = unum_node_pre;
    elseif flag == 2029
        coeff_store{3,1} = flag;
        coeff_store{3,2} = unum_node_pre;
    elseif flag == 2648
        coeff_store{4,1} = flag;
        coeff_store{4,2} = unum_node_pre;
    end

    %% Output real-time monitoring data
    fprintf('T_now = %.10f, Linf_now = %.8f, dt_now = %.8g, mu_now = %.8f, flag = %d\n', ...
        TT(flag), U_max(flag), dt_now, Mu(flag), flag-1)
end

%% Save the numerical results
filename = ['BioMT_TSM4_s=' num2str(s) '_p=' num2str(p) '_N=' num2str(N) '_Bound=' num2str(U_max_bound) '.mat'];
save(filename,"s","p","gamma", ...
    "N","mu_initial","dt_initial", ...
    "U_max_bound","eps_time","safe_factor","eps_space","ratio", ...
    "Mu","TT","Findex","U_max","Ugrad_L2","Mass","Energy","time_iter","coeff_store");