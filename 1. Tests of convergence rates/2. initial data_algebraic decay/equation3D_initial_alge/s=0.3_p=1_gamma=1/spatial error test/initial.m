% This function is used in alge3D_error_space.m to represent the initial
% value function \psi_0(x,y,z).

function f = initial(x,y,z)
f = 1./(1+x.^2+y.^2+z.^2).^(1+sqrt(2));