% This function is used in generate_ref_1D_alge.m to represent the initial
% value function \psi_0(x).

function f = initial(x)
f = 1./(1+x.^2).^sqrt(2);