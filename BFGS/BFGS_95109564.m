function [xopt, fopt, iter] = BFGS_95109564 (f, gf, x0, Stop_tol, fbar, c1, c2, varargin)
ff  = @(x) f(x, varargin{:});
gff = @(x) gf(x, varargin{:});
[xopt, fopt, iter] = mainBFGS (ff, gff, x0, Stop_tol, fbar, c1, c2);
end

function [xopt, fopt, iter] = mainBFGS (ff, gff, x0, Stop_tol, fbar, c1, c2)
addpath('./counters')
tol = inf;
x = x0;
iter = 0;
alpha = 1;
n = length(x0);
C = eye(n, n);
I = eye(n, n);
start_flag = 1;
gk = inf;
gk = gff(x);
while 1
    norm(gk)
    if (norm(gk) < Stop_tol)
        break;
    end
    
    p = -C*gk;
    LineFunction = @(a) ff(x + a.*p);
    
    LineFunctionDerv = @(a) gff(x + a.*p)'*p;
    
    alpha = bracketing_95109564(LineFunction , LineFunctionDerv, 1, fbar, c1, c2);
    
    x = x + alpha.*p;
    gkpast = gk;
    gk = gff(x);
    delta = alpha.*p;
    gamma = gk -gkpast;
    a = 1/(delta'*gamma);
    if start_flag == 1
        C = (gamma'*delta)/(gamma'*gamma)*I;
        start_flag = 0;
    else
        C = (I - a *delta*gamma')*C*(I - a * gamma * delta') + a * delta *delta';
    end
    
    iter = iter + 1;
    %ff(x)
end
xopt = x;
fopt = ff(xopt);
end
