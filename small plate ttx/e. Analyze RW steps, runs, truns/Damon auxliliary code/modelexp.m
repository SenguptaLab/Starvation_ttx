function y=modelexp(p,x);

% models an exponential decay

y=p(1)*exp(-x/p(2));