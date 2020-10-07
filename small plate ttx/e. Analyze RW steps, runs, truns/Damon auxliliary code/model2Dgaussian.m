function y=model2Dgaussian(p,x);

% makes a gaussian.  real(x) are x values, imag(x) are y values
% p = [baseline height xcenter ycenter xsigma ysigma rotation]

x0=real(x);
y0=imag(x);

x1=(x0-p(3))*cos(p(7))-(y0-p(4))*sin(p(7));   % rotate it
y1=(x0-p(3))*sin(p(7))+(y0-p(4))*cos(p(7));

y=p(1)+p(2)*exp(- x1.^2/p(5)^2 - y1.^2/p(6)^2);