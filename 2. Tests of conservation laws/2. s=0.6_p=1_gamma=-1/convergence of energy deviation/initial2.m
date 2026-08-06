% This function is used in MT_TSM4_energy_order2D.m to represent the
% initial value function \psi_0(x,y).

function f = initial2(x,y)
f = 1./(1+x.^2+y.^2).^(1/2+sqrt(2));