% This function is used in gauss3D_error_time.m to represent the initial
% value function \psi_0(x,y,z).

function f = initial(x,y,z)
f = exp(-x.^2-y.^2-z.^2);