% This function is used in equation2D_collision.m to represent the initial
% value function \psi_0(x,y).

function f = initial(x,y)
f = exp(-(x+3).^2-y.^2).*exp(1i*5*x) + exp(-(x-3).^2-y.^2).*exp(-1i*5*x);