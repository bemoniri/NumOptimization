function [alpha_star] = bracketing_95109564(phi, phiprime, alpha1, fbar, c1, c2)

alpha(1) = 0;   %previous alpha!
alpha(2) = alpha1; % current alpha!

phival = phi(alpha(2)); % current phi
lastphival = phi(alpha(1)); % last phi

phizero = phi(0);
phiprimezero = phiprime(0);
alpha_max = (fbar - phizero)/(c1*phiprimezero);

i = 1;
while true    
    if(phival <= fbar)
        alpha_star = alpha(2);        
        break
    end
    
    if( (phival > phizero + c1*alpha(2)*phiprimezero) ||( (phival > lastphival) && (i > 1 )) )
        alpha_star = sectioning_95109564(phi, phiprime, c1, c2, alpha(1), alpha(2));
        break
    end
    
    phiprimeval = phiprime(alpha(2));
    
    if ( abs(phiprimeval) <=  -c2*phiprimezero )
        alpha_star = alpha(2);
        break
    end
    
    if ( phiprimeval  >= 0)
        alpha_star = sectioning_95109564(phi, phiprime, c1, c2, alpha(2), alpha(1));
        break
    end
    
    if(alpha_max <= 2*alpha(2) - alpha(1))
        alpha(1) = alpha(2);
        alpha(2) = alpha_max;
    else
        A = 2*alpha(2) - alpha(1);
        B = min(alpha_max, alpha(2) + 10*(alpha(2) - alpha(1)));
        alpha(1) = alpha(2);
        alpha(2) = 0.5*(A+B);
    end
        
    i = i + 1;
    %if i == maxit
    %    disp('Maximum number of iteration for Line Search reached');
    %    alpha_star = alpha;
    %    return;
    %end
end
end