% This function is used in equation3D_diffusion.m to represent the initial
% value function \psi_0(x,y,z).

function f = initial(x,y,z)
f = exp(-4*(x-1).^2 - 4*(y-0.25).^2 - 4*(z-2.5).^2) ...
    + exp(-(x+2.5).^2 - (y+0.25).^2 - (z+1.5).^2) ...
    + exp(-2*(x-3).^2 - 2*(y+0.5).^2 - 2*(z+1.5).^2);