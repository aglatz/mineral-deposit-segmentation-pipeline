function [y, x] = chiqqplot_mod(y,df,class)
% Modified version of 'chiplot()' from the LIBRA library

%CHIQQPLOT produces a Quantile-Quantile-plot of the vector y 
% versus the square root of the quantiles of the chi-squared distribution.
%
% Required input arguments:
%       y  : row or column vector 
%       p  : degrees of freedom of the chi-squared distribution 
%
% Optional input argument:
%    class : a string used for the y-label and the title(default: ' ')
%
% I/O: chiqqplot(y,p,class)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
%Written by Nele Smets
% Last update: 23/10/2003

set(gcf,'Name', 'Chisquare QQ-plot', 'NumberTitle', 'off')
n=length(y);
if nargin==2
    class='';
end
scale = iqr(y).^2 / (chi2inv(.75, df) - chi2inv(.25, df));
for i=1:n
    p = (i - .5)/n; %(i-1/3) / (n+1/3);
	x(i) = chi2inv(p, df);
    g = (x(i)^(-1+df/2)) / (exp(x(i)/2)*gamma(df/2)*(sqrt(2)^df));
    se(i) = (scale/g)*sqrt(p*(1-p)/n);
end
x=sqrt(x);
se=sqrt(se);
y=sort(y);
plot(x, x, 'k');
hold on;
plot(x, x - 1.95*se, '--k', x, x + 1.95*se, '--k');
plot(x, y, 'm');
xlabel('\bf Square root of the quantiles of the chi-squared distribution');
% if strcmp(class,'MCDCOV')
	ylabel('\bf Square root of the robust distances');
% elseif strcmp(class,'COV')
%    ylabel('\bf Mahalanobis distance');
% else
%     ylabel('Distance');
% end
% title(class);
