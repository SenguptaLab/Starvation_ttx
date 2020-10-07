function [A, lambda] = my_exp_fit4(x,y)

indx = find(y>0);
fx = x(indx);
fy = log(y(indx)); 
[p,s] = polyfit(fx,fy,1); % fit to line
A = exp(p(2));
lambda = -1/p(1); 

return;