function [result, N] = GSS (f, a, b, tol, varargin)
    N = ceil(log(tol/2)./log(1 - (3-sqrt(5))/2));    
    ff = @(x) f(x, varargin{:});    
    [result] = mainGSS (ff, a, b, tol, 0, 0);
end

function [result] = mainGSS (f, a, b, tol, numbwe, value)    
    
    rho = (3-sqrt(5))./2;
    h = b - a;
    
    a0 = a;
    b0 = b;
    a1 = a0 + rho.*h;
    b1 = a0 + (1-rho).*h;
    
    if (numbwe == 0)
        fb1 = f(b1);
        fa1 = f(a1);
    end
    
    if (numbwe == 1)
        fb1 = value;
        fa1 = f(a1);
    end
    
    if(numbwe == 2)
        fb1 = f(b1);
        fa1 = value;
    end
    
    
    if (b-a < tol)
       result = (b+a)/2;
       return
    end
    
    if (fa1<fb1)
       result = mainGSS(f, a0, b1, tol, 1, fa1);
    end
    
    if (fa1>=fb1)
       result = mainGSS(f, a1, b0, tol, 2, fb1);      
    end
end
