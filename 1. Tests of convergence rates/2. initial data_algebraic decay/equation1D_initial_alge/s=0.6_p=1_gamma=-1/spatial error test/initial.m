% This function is used in alge1D_error_space.m to represent the initial
% value function \psi_0(x).

function f = initial(x)
f = 1./(1+x.^2).^sqrt(2);