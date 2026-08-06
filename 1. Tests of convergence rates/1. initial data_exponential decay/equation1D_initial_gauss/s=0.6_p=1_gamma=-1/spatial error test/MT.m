% This program returns the value of the k-th MT function at point x, where
% MT(k,x) = (1i)^k*sqrt(2/pi)*(1+2*i*x)^k/(1-2*i*x)^(k+1).

function mt = MT(k,x)
mt = (1i).^k.*sqrt(2/pi).*((1+2*1i*x)./(1-2*1i*x)).^k./(1-2*1i*x);