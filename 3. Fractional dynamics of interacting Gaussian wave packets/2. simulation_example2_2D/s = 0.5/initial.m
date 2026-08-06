% This function is used in equation2D_diffusion.m to represent the initial
% value function \psi_0(x,y).

function f = initial(x,y)
f = exp(-8*(x-3).^2-8*(y-2).^2) + exp(-(x+4).^2-(y+4).^2) + exp(-2.5*(x-3).^2-2.5*(y+4).^2);