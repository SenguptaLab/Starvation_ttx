function y=modelexp_norm(p,x);

% models an exponential decay

y=1/p(1)*exp(-x/p(1));