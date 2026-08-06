% This function is used in equation2D_dynamic.m to represent the initial
% value function \psi_0(x,y).

function f = initial(x,y)
% f = exp(-1/2*x.^2-1/2*y.^2);
% f = exp(-x.^2-y.^2);
% f = 1/2*exp(-1/2*x.^2-1/5*y.^2);
% f = 2./(1+x.^2+y.^2).^3;
% f = exp(-x.^2-y.^2).*exp(1i*5*x+1i*3*y);
% f = exp(-x.^2-y.^2).*exp(2*1i*(x.^2+y.^2));
% f = exp(-(x-2).^2-y.^2) - exp(-(x+2).^2-y.^2);
% f = exp(-(x-3).^2-2*y.^2).*exp(-1i*5*x) + exp(-(x+3).^2-2*y.^2).*exp(1i*5*x);
f = exp(-(x+3).^2-3*(y-0.5).^2).*exp(1i*4*x) + exp(-(x-3).^2-3*(y+0.5).^2).*exp(-1i*4*x);
% f = exp(-(x-2).^2-y.^2) + exp(-(x+2).^2-y.^2)*exp(1i*pi*3/2);
% f = exp(-(x-5).^2-(y-5).^2) + exp(-(1/2*x+3).^2-(1/2*y+3).^2) + exp(-(1/3*x-2).^2-(1/3*y+4).^2);
% f = exp(-(1/2*x-4).^2-(y-4).^2) + exp(-(2*x+2).^2-(1/3*y+2).^2);
% f = 1./(1+x.^2+y.^2).^(1/2+sqrt(2));
% f = exp(-x.^2-y.^2).*sin(x+x.*y+y.^2);
% f = exp(-(x-5).^2-(y-5).^2) + exp(-(1/2*x+3).^2-(1/2*y+3).^2);
% f = exp(-(1/2*x-5/2).^2-(y-3/2).^2) + exp(-(2*x+2).^2-(1/3*y+2).^2);
% f = 1./(1+x.^2+y.^2).^(1+sqrt(2)/10);