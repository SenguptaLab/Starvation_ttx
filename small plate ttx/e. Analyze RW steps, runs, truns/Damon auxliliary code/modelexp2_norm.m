function y=modelexp2_norm(p,x);

% models an exponential decay

y=p(3)/p(1)*exp(-x/p(1))+(1-p(3))/p(2)*exp(-x/p(2));