function y=modelstep(p,x);

% models the step response as the sum of two exponentials
c0=find(x<0);
c1=find(x>=0);

y(c0)=p(1);
if 1
    y=p(1)+p(2).*(exp(-x(c1)./p(3)-exp(-x(c1)./p(4));
else
    y=p(1)+p(2).*(1-exp(-x(c1)/p(3))).*exp(-x(c1)/p(4));
end