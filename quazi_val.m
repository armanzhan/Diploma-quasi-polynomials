function [y] = quazi_val(p,q,tau,x)
y = polyval(p,x) + polyval(q,x) .* exp(tau * x);
end