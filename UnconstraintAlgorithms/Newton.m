function [xopt, fopt, iter] = Newton (f, gf, Hf, x0, Stop_tol, GSS_tol, varargin)    
    ff  = @(x) f(x, varargin{:});    
    gff = @(x) gf(x, varargin{:});
    Hff = @(x) Hf(x, varargin{:});
    [xopt, fopt, iter] = mainNewton (ff, gff,Hff, x0, Stop_tol, GSS_tol);
end

function [xopt, fopt, iter] = mainNewton (ff, gff, Hff, x0, Stop_tol, GSS_tol)        
    
    tol = inf;
    x = x0;
    iter = 0;
    while (tol>Stop_tol)        
        Hessian = Hff(x);
        %E = eig(Hessian);
        %if(E(1) < 0)
        %    Hessian = Hessian - (-10+E(1))*eye(size(Hessian));
        %    disp('besmel')
        %end
        %eig(Hessian)
        p = -linsolve(Hessian,gff(x)')'; %Direction of Movement            
        LineFunction = @(a) ff(x + a.*p);
        alpha = GSS(LineFunction, 0, 100, GSS_tol);
        x = x + alpha.*p;        
        tol = norm(p,2);
        iter = iter + 1;
    end
    xopt = x;
    fopt = ff(xopt);
end
