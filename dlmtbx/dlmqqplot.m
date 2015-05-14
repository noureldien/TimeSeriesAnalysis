function out=dlmqqplot(dlm,ind)
%DLMQQPLOT normal probability plot for DLM fit residuals 

if nargin<2, ind = 1; end
% p = size(dlm.F,1); % number of series

y = dlm.resid2(:,ind);
igood = not(isnan(y));
yy = y(igood);
h=qqplot(yy);
%set(h,'markerfacecolor','black');
title('Normal probability plot for the residuals')
grid on;

if nargout>0
  out=h;
end
