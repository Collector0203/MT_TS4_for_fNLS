% This function is used in generate_ref_2D_alge.m to represent the initial
% value function \psi_0(x,y).

function f = initial(x,y)
f = 1./(1+x.^2+y.^2).^(1/2+sqrt(2));