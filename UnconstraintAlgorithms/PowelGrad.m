function [out] = PowelGrad(x)
out = [2*(x(1)+10*x(2)) + 40*(x(1)-x(4))^3; ...
            20*(x(1)+10*x(2)) + 4*(x(2)-2*x(3))^3; ...
            10*(x(3)-x(4)) - 8*(x(2)-2*x(3))^3; ...
            -10*(x(3)-x(4)) - 40*(x(1)-x(4))^3];

incgradcount;
end

