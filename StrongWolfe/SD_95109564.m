function [xopt, fopt, iter] = SD_95109564 (f, gf, x0, Stop_tol, fbar, c1, c2, varargin)
ff  = @(x) f(x, varargin{:});
gff = @(x) gf(x, varargin{:});
[xopt, fopt, iter] = mainSD (ff, gff, x0, Stop_tol, fbar, c1, c2);
end

function [xopt, fopt, iter] = mainSD (ff, gff, x0, Stop_tol, fbar, c1, c2)
addpath('./counters')
tol = inf;
x = x0;
iter = 0;
alpha = 1;

while (tol>Stop_tol)
    
    p = -gff(x);
    LineFunction = @(a) ff(x + a.*p');
    LineFunctionDerv = @(a) gff(x + a.*p')'*p;   
    
    %alpha = GSS(LineFunction, 0, 100, GSS_tol);
    f_bar = 0.1*LineFunction(0);
    alpha = bracketing_95109564(LineFunction , LineFunctionDerv, 1, fbar, c1, c2);                
    %alpha = lsa(ff, gff,x,d,a1);
    
    %if (~checkstrongwolfe(LineFunction, LineFunctionDerv, alpha, c1, c2))
    %    disp('NOT STRONG WOLFE')
    %    %break
    %end
    
    x = x + alpha.*p';
    tol = norm(alpha.*p,2);    
    iter = iter + 1;
   ff(x)
end
xopt = x;
fopt = ff(xopt);
end
