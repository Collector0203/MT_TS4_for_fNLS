% This function returns the exact solution of 1i*\partial_t \psi(x,t) =
% 2*gamma*|\psi(x,t)|^{2*p}*\psi(x,t), where x \in R^d, t \in R. 
% 
% input: 
% 1. the values of \psi(x,0) "psi".
% 2. time step "dt".
% 3. nonlinearity parameters "p" and "gamma".
% 
% output: the values of \psi(x,dt).

function Psi = subNonlinear(psi,dt,gamma,p)
Psi = exp(-2*1i*gamma*dt*abs(psi).^(2*p)).*psi;