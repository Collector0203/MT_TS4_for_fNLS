% This function is used in gauss1D_error_space.m to represent the initial
% value function \psi_0(x).

function f = initial(x)
f = exp(-x.^2);