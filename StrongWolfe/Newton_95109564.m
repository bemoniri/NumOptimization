function [xopt, fopt, iter] = Newton_95109564 (f, gf, Hf, x0, Stop_tol, fbar, c1, c2, varargin)
ff  = @(x) f(x, varargin{:});
gff = @(x) gf(x, varargin{:});
Hff = @(x) Hf(x, varargin{:});
[xopt, fopt, iter] = mainNewton (ff, gff,Hff, x0, Stop_tol, fbar, c1, c2);
end

function [xopt, fopt, iter] = mainNewton (ff, gff, Hff, x0, Stop_tol, fbar, c1, c2)
addpath('./counters')
tol = inf;
x = x0;
iter = 0;
while (tol>Stop_tol)
    Hessian = Hff(x);
    p = -linsolve(Hessian,gff(x))'; %Direction of Movement
    LineFunction = @(a) ff(x + a.*p);
    LineFunctionDerv = @(a) p*gff(x + a.*p);        
    alpha = bracketing_95109564(LineFunction , LineFunctionDerv, 1, fbar, c1, c2);        
    x = x + alpha.*p;
    tol = norm(p,2);
    iter = iter + 1;
    %ff(x)
end
xopt = x;
fopt = ff(xopt);
end
