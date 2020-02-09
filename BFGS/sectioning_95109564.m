function output = sectioning_95109564 (phi, phiprime, c1, c2, a, b)

phi0 = phi(0);
phiprime0 = phiprime(0);
i = 1;
alpha = 1;
alpha
while true
    phi_a = phi(a);
    phiprime_a = phiprime(a);
    temp_phiprime_a = (b-a)*phiprime_a;
    
    zeta = -temp_phiprime_a/(2*phi_a);
    phi_zeta = phi(a + zeta*(b-a));       
    
    if(phi_zeta < phi_a)
        alpha = a + zeta.*(b-a);
    else
        alpha = 0.5*(a+b);
    end
    
    %alpha = 0.5*(a+b);   % Choose alpha based on a and b    
    phival = phi(alpha); % Evaluate phi(alpha)
    
    if ( (phival > phi0 + c1*alpha*phiprime0) || (phival >= phi_a))
        b = alpha;        
    else
        phiprimeval = phiprime(alpha); % Evaluate phi'(alpha)
        
        if (abs(phiprimeval) <= -c2*phiprime0)
            output = alpha;
            break
        end
        
        if ((b-a)*phiprimeval >= 0)
            b = a;
        end
        a = alpha;
    end   
    i = i +1;
    if (i > 10000)
        output = alpha;
        break
    end
end
end