function y=model2Dgaussian_double(p,x);

% makes a gaussian.  real(x) are x values, imag(x) are y values
% p = [baseline0 baseline1 height0 height1 xcenter0 ycenter0 xcenter1 ycenter1 xsigma ysigma rotation]
%%%%      1         2         3       4       5        6        7        8
%%%%      9         10        11

x01=real(x(1:length(x)/2));
y01=imag(x(1:length(x)/2));
x02=real(x(length(x)/2+1:end));
y02=imag(x(length(x)/2+1:end));

x11=(x01-p(5))*cos(p(11))-(y01-p(6))*sin(p(11));   % rotate it
y11=(x01-p(5))*sin(p(11))+(y01-p(6))*cos(p(11));
x12=(x02-p(7))*cos(p(11))-(y02-p(8))*sin(p(11));   % rotate it
y12=(x02-p(7))*sin(p(11))+(y02-p(8))*cos(p(11));

y(1:length(x)/2)=p(1)+p(3)*exp(- x11.^2/p(9)^2 - y11.^2/p(10)^2);
y(length(x)/2+1:length(x))=p(2)+p(4)*exp(- x12.^2/p(9)^2 - y12.^2/p(10)^2);
y=y';