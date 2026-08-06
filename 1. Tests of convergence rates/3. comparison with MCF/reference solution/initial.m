% This function is used in generate_ref_comp_MCF.m to represent the initial
% value function \psi_0(x,y).

function f = initial(x,y)
f = sech(x).*sech(y).*exp(1i*(x+y));