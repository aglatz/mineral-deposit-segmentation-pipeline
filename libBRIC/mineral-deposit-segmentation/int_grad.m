function [res] = int_grad(x, dx, x0, gz, t)
M = gz~=0;
P = polyfit(x(M), gz(M), 3);
x_new = min(x)-dx/2:0.1:max(x)+dx/2;
gz_new = polyval(P, x_new);
x_norm = (x_new-x0)/(dx/2);
sfc = sinc(x_norm);
sfc = sfc./((1-x_norm.^2));
sfc(isinf(sfc)) = 0.5;
% sfc = double(abs(x_norm) < 1); 
res = trapz(x_new, exp(sqrt(-1)*gz_new*t).*sfc);