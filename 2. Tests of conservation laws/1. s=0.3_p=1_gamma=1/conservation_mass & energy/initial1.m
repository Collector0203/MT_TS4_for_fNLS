% This function is used in MT_TSM4_conservation2D.m to represent the
% initial value function \psi_0(x,y).

function f = initial1(x,y)
f = exp(-x.^2-y.^2);