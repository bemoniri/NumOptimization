function [xopt, fopt, iter] = SD (f, gf, x0, Stop_tol, GSS_tol, varargin)    
    ff  = @(x) f(x, varargin{:});    
    gff = @(x) gf(x, varargin{:});
    [xopt, fopt, iter] = mainSD (ff, gff, x0, Stop_tol, GSS_tol);
end

function [xopt, fopt, iter] = mainSD (ff, gff, x0, Stop_tol, GSS_tol)            
    tol = inf;
    x = x0;
    iter = 0;    
    while (tol>Stop_tol)
        p = -gff(x);     
        LineFunction = @(a) ff(x + a.*p);
        alpha = GSS(LineFunction, 0, 100, GSS_tol);
        x = x + alpha.*p;
        tol = norm(alpha.*p,2);
        iter = iter + 1;      
    end
    xopt = x;
    fopt = ff(xopt);
end
