% This function is used in generate_ref_3D_gauss.m to represent the initial
% value function \psi_0(x,y,z).

function f = initial(x,y,z)
f = exp(-x.^2-y.^2-z.^2);