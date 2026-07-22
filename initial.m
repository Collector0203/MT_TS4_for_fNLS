% This function is used in generate_ref_1D_gauss.m to represent the initial
% value function \psi_0(x).

function f = initial(x)
f = exp(-x.^2);