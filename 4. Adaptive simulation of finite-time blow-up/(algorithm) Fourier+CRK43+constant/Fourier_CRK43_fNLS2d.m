%% Fourier_CRK43_fNLS2d.m
% Fourier pseudospectral + Driscoll CRK43 solver for the 2D fractional NLS
%
%   i * psi_t = (1/2) * (-Delta)^s psi + gamma * |psi|^(2p) * psi
%
% on the periodic square [-D,D]^2.
%
% Required user-supplied file:
%   initial.m
% with the interface
%   f = initial(X,Y)
% where X and Y are N-by-N grids and f is the initial value on that grid.
%
% Time discretization:
%   Composite Runge--Kutta method CRK43 of Driscoll (2002).
%   Low Fourier modes are advanced by classical explicit RK4.
%   High Fourier modes have their linear part advanced by the linearly
%   implicit RK formula used in CRK43, while their nonlinear part remains
%   explicit.  If there are no high modes, the method reduces exactly to
%   classical RK4.
%
% Notes:
%   1. The Fourier multiplier for (-Delta)^s in 2D is
%          (kx^2 + ky^2)^s.
%   2. Since the physical domain is [-D,D]^2, the Fourier wavenumbers are
%          (pi/D) * integer_modes.
%   3. This code uses a fixed time step dt throughout the computation.
%   4. The nonlinear term is evaluated pseudospectrally without dealiasing.

clear; clc;

%% -------------------- Equation parameters --------------------

s     = 0.7;
p     = 1.5;
gamma = -1;                  % focusing: gamma = -1; defocusing: gamma = +1
T     = 1.3;                % terminal time

%% -------------------- Discretization parameters --------------------

N = 1024;
Nreal  = 2*N;                % Nreal Fourier modes in x and N Fourier modes in y
D  = 10*pi;                  % physical domain is [-D,D]^2
dt = 1e-2;                   % fixed time step

%% -------------------- Parameter checks --------------------

if ~isscalar(s) || ~isreal(s) || s <= 0
    error('s must be a positive real scalar.');
end
if ~isscalar(p) || ~isreal(p) || p <= 0
    error('p must be a positive real scalar.');
end
if ~isscalar(gamma) || ~isreal(gamma)
    error('gamma must be a real scalar.');
end
if ~isscalar(T) || ~isreal(T) || T <= 0
    error('T must be a positive real scalar.');
end
if ~isscalar(D) || ~isreal(D) || D <= 0
    error('D must be a positive real scalar.');
end
if ~isscalar(dt) || ~isreal(dt) || dt <= 0
    error('dt must be a positive real scalar.');
end
if ~isscalar(Nreal) || ~isreal(Nreal) || Nreal <= 0 || Nreal ~= round(Nreal) || mod(Nreal,2) ~= 0
    error('N must be a positive even integer.');
end

Nt_real = T / dt;
Nt = round(Nt_real);
if abs(Nt_real - Nt) > 100 * eps(max(1, abs(Nt_real)))
    error('T/dt must be an integer to numerical precision.');
end
TT = (0:Nt)*dt;              % discrete time

%% -------------------- Periodic Fourier grid on [-D,D]^2 --------------------

dx = 2*D / Nreal;
dy = dx;
dxdy = dx * dy;

% Endpoint excluded periodic grid.
x = -D + dx*(0:Nreal-1);
y = -D + dy*(0:Nreal-1);

% First array index represents x, second array index represents y.
[X,Y] = ndgrid(x,y);

%% -------------------- Initial data supplied by initial.m --------------------

if exist('initial', 'file') ~= 2
    error('initial.m was not found on the current MATLAB path.');
end

psi = initial(X,Y);

if ~isequal(size(psi), [Nreal,Nreal])
    error('initial.m must return an N-by-N array f.');
end
if ~isnumeric(psi) || any(~isfinite(psi(:)))
    error('The initial value returned by initial.m must be finite and numeric.');
end

psi = double(psi);
clear X Y

%% -------------------- Fourier wavenumbers and linear symbol --------------------

% Fourier frequencies on an interval of length 2D.
mode = [0:Nreal/2-1, -Nreal/2:-1];
k = (pi/D) * mode;

% The arrays below have size N-by-N and correspond to ndgrid ordering.
KX = k(:);
KY = k(:).';

% Fourier multiplier of (-Delta)^s.
kpow = (KX.^2 + KY.^2).^s;

% Linear multiplier in psi_t = L psi + nonlinear term.
Lhat = -0.5i * kpow;

%% -------------------- Driscoll CRK43 mode partition --------------------

% Driscoll's slow-mode criterion is |dt * Lhat| < 2.8.
rk4ImagCutoff = 2.8;
rho = abs(dt * Lhat);
lowMask = (rho < rk4ImagCutoff);
highMask = ~lowMask;

nLow = nnz(lowMask);
nHigh = nnz(highMask);

% Linear operator acting on high modes only.
Lambda = Lhat .* highMask;

% Diagonal implicit denominators in the CRK43 high-mode stages.
den2 = 1 - (dt/3) * Lambda;
den3 = 1 - dt * Lambda;

fprintf('Fourier_CRK43_fNLS2d: N = %d, total modes = %d, Nt = %d\n', Nreal, Nreal^2, Nt);
fprintf('Domain: [-%.8f, %.8f]^2, fixed dt = %.8e\n', D, D, dt);
fprintf('CRK43 split: low modes = %d, high modes = %d\n', nLow, nHigh);
if nHigh == 0
    fprintf('All modes are low modes: CRK43 reduces exactly to classical RK4.\n');
else
    fprintf('Full CRK43 is active: high-mode linear terms use the implicit stages.\n');
end

%% -------------------- Initial Fourier representation --------------------

Uhat = fft2(psi);

%% -------------------- Diagnostic storage --------------------

Mass         = zeros(Nt+1,1);
Energy       = zeros(Nt+1,1);
relMassErr   = zeros(Nt+1,1);
relEnergyErr = zeros(Nt+1,1);
U_max        = zeros(Nt+1,1);

% Initial diagnostics at t = 0.
absPsi = abs(psi);
U_max(1) = max(absPsi(:));
Mass(1) = dxdy * sum(absPsi(:).^2);

kinetic = 0.5 * dxdy/(Nreal^2) * sum(kpow(:) .* abs(Uhat(:)).^2);
potential = gamma/(p+1) * dxdy * sum(absPsi(:).^(2*p+2));
Energy(1) = kinetic + potential;

M0 = Mass(1);
E0 = Energy(1);

if M0 == 0
    error('The initial discrete mass is zero; relMassErr is undefined.');
end
if E0 == 0
    warning(['The initial discrete energy is exactly zero. ', ...
             'relEnergyErr is undefined and will be stored as NaN after t=0.']);
    energyRelativeDefined = false;
    relEnergyErr(1) = 0;
else
    energyRelativeDefined = true;
    if abs(E0) <= 100*eps(max(1,abs(E0)))
        warning('The initial discrete energy is very small; relEnergyErr may be ill-conditioned.');
    end
end

fprintf('t = %.8f, U_max = %.8f\n', 0, U_max(1));

%% -------------------- Fixed-step Driscoll CRK43 time integration --------------------

tic

for n = 1:Nt

    % Partition the current Fourier coefficients into high and low modes.
    Yh = highMask .* Uhat;
    Zl = lowMask  .* Uhat;

    % ================================================================
    % Stage 1
    % Y1 = Yh, Z1 = Zl.  The physical-space value psi is already the
    % inverse FFT of the current Uhat from the previous diagnostic step.
    % ================================================================
    Nphys = -1i * gamma * abs(psi).^(2*p) .* psi;
    Nhat = fft2(Nphys);

    f1 = highMask .* Nhat;
    g1 = lowMask .* (Lhat .* Zl + Nhat);

    % ================================================================
    % Stage 2
    % ================================================================
    Y2 = highMask .* ((Yh + 0.5*dt*f1 + (dt/6)*Lambda.*Yh) ./ den2);
    Z2 = lowMask  .* ( Zl + 0.5*dt*g1 );

    % Accumulators for the final CRK43 weighted combination.
    Yacc = f1 + Lambda.*Yh;
    Zacc = g1;

    psiStage = ifft2(Y2 + Z2);
    Nphys = -1i * gamma * abs(psiStage).^(2*p) .* psiStage;
    Nhat = fft2(Nphys);

    f2 = highMask .* Nhat;
    g2 = lowMask .* (Lhat .* Z2 + Nhat);

    % ================================================================
    % Stage 3
    % ================================================================
    Y3 = highMask .* ((Yh + 0.5*dt*f2 + 0.5*dt*Lambda.*Yh ...
                      - dt*Lambda.*Y2) ./ den3);
    Z3 = lowMask  .* ( Zl + 0.5*dt*g2 );

    Yacc = Yacc + 2*(f2 + Lambda.*Y2);
    Zacc = Zacc + 2*g2;

    clear f1 g1 f2 g2 Y2 Z2

    psiStage = ifft2(Y3 + Z3);
    Nphys = -1i * gamma * abs(psiStage).^(2*p) .* psiStage;
    Nhat = fft2(Nphys);

    f3 = highMask .* Nhat;
    g3 = lowMask .* (Lhat .* Z3 + Nhat);

    % ================================================================
    % Stage 4
    % ================================================================
    Y4 = highMask .* ((Yh + dt*f3 + (2*dt/3)*Lambda.*Y3) ./ den2);
    Z4 = lowMask  .* ( Zl + dt*g3 );

    Yacc = Yacc + 2*(f3 + Lambda.*Y3);
    Zacc = Zacc + 2*g3;

    clear f3 g3 Y3 Z3

    psiStage = ifft2(Y4 + Z4);
    Nphys = -1i * gamma * abs(psiStage).^(2*p) .* psiStage;
    Nhat = fft2(Nphys);

    f4 = highMask .* Nhat;
    g4 = lowMask .* (Lhat .* Z4 + Nhat);

    % ================================================================
    % Final CRK43 update
    % ================================================================
    Yacc = Yacc + f4 + Lambda.*Y4;
    Zacc = Zacc + g4;

    Uhat = Uhat + (dt/6) * (Yacc + Zacc);

    clear Yh Zl Y4 Z4 Yacc Zacc f4 g4 psiStage Nphys Nhat

    %% ---------------- Diagnostics at t = n*dt ----------------

    psi = ifft2(Uhat);
    absPsi = abs(psi);

    U_max(n+1) = max(absPsi(:));
    Mass(n+1) = dxdy * sum(absPsi(:).^2);

    kinetic = 0.5 * dxdy/(Nreal^2) * sum(kpow(:) .* abs(Uhat(:)).^2);
    potential = gamma/(p+1) * dxdy * sum(absPsi(:).^(2*p+2));
    Energy(n+1) = kinetic + potential;

    relMassErr(n+1) = abs(Mass(n+1) - M0) / abs(M0);
    if energyRelativeDefined
        relEnergyErr(n+1) = abs(Energy(n+1) - E0) / abs(E0);
    else
        relEnergyErr(n+1) = NaN;
    end

    % Required real-time output at every discrete time level.
    fprintf('t = %.8f, U_max = %.8f\n', n*dt, U_max(n+1));
end

toc

%% -------------------- Save required output variables --------------------

% % Output file.  The required variables are saved to this MAT-file.
% outputFile = ['Fourier_CRK43_s=',num2str(s),'_p=',num2str(p),'_N=',num2str(N),'.mat'];
% 
% save(outputFile, ...
%     's', 'p', 'gamma', 'T', 'N', 'dt', 'D', 'TT', ...
%     'Mass', 'Energy', 'relMassErr', 'relEnergyErr', 'U_max');
% 
% fprintf('Saved results to %s\n', outputFile);