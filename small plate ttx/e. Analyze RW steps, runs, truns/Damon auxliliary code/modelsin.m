function y=modelsin(p,x);

% model of sine, expected in various angular fits i'm doing

y=p(1)+p(2)*sin(x+p(3));