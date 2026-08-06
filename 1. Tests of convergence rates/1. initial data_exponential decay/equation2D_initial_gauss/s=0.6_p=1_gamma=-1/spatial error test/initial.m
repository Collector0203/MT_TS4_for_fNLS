% This function is used in gauss2D_error_space.m to represent the initial
% value function \psi_0(x,y).

function f = initial(x,y)
f = exp(-x.^2-y.^2);