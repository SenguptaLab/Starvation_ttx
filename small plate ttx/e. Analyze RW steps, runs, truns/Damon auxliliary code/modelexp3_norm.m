function y=modelexp3_norm(p,x);

% models an exponential decay

y=p(4)/p(1)*exp(-x/p(1))+p(5)/p(2)*exp(-x/p(2))+(1-p(4)-p(5))/p(3)*exp(-x/p(3));