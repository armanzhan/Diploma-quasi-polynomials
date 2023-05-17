function [y] = quazi_val(p,q,tau,x) % аналог функции polyval для квазиполиномов с одним запаздыванием
y = polyval(p,x) + polyval(q,x) .* exp(tau * x);
end
