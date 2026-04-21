function [an] = integral_an(v,n)
fun = @(x) (exp(-(1-cos(pi*x))/(2*pi*v))).*cos(n*pi*x);
an=2*integral(fun,0,1);
end