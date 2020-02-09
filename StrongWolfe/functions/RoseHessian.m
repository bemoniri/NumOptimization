function [Hessian] = RoseHessian(x)
    Hessian = [1200.*x(1)-400.*x(2)+2, -400.*x(1); -400.*x(1), 200];
    inchesscount;
end

