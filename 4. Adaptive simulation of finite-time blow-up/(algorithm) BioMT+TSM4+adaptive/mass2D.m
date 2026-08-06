% This function computes the discrete mass of u(x,y) by the values of u on
% the MT-nodes.
% 
% input: 
% 1. the function value matrix of u(x,y) on MT-nodes "U".
% 2. truncation parameter "N".
% 3. scaling factor "mu".
% 4. common matrix "opt_CosTheta", which is determined by truncation
% parameter "N".
% 
% output: discrete mass of u(x,y).
%
% (This function has been optimized for efficiency.)

function mass = mass2D(U,N,mu,opt_CosTheta)
mass = (pi/4/mu/N)^2*sum(abs(U./opt_CosTheta).^2,'all');