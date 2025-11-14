function [a0] = integral_a0(v)
fun = @(x) exp(-(1-cos(pi*x))/(2*pi*v));
a0=integral(fun,0,1);
end