function [temp] = checkstrongwolfe(phi, phiprime, alpha, c1, c2)
temp = 0;
if (phi(alpha) <= phi(0) + c1.*alpha*phiprime(0)) % Armijo
   if (abs(phiprime(alpha)) <= c2.*abs(phiprime(0))) % Second Strong Wolfe
      temp = 1; 
   end    
end
end

