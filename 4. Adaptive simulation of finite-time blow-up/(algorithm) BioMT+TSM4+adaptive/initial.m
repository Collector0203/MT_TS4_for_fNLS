% This function is used in equation2D_blow_gauss.m to represent the initial
% value function \psi_0(x,y).

function f = initial(x,y)
f = exp(-1/2*x.^2-1/2*y.^2);