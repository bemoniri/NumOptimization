function [x_min, N] = GSS(f, a, b, epsilon, varargin)

rho = (3-sqrt(5))/2;
N = 0;

while b-a > epsilon
    N = N+1;
    x = (1-rho) *a + b*rho;
    y = (1-rho) *b + a*rho;
    if f(x, varargin{:})<f(y, varargin{:})
        b = y;
    else
        a = x;
    end
end

x_min = (a+b)/2;
