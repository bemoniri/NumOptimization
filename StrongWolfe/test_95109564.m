clc
clear all
format long
addpath('./functions')

%% SD - Powel Function

global gradcount;
gradcount = 0;
global funccount;
funccount = 0;
global hesscount;
hesscount = 0;

f = @(x) Powel(x);
gf = @(x) PowelGrad(x);

fbar = 0.01;
c1 = 0.0001;
c2 = 0.1;
stop_tol = 1e-3;

[xopt, fopt, iter] = SD_95109564(f, gf, [1,2,2,2], stop_tol,  fbar, c1 ,c2)

gradcount
funccount
hesscount

%% Newton - Powel Function
clc
global gradcount;
gradcount = 0;
global funccount;
funccount = 0;
global hesscount;
hesscount = 0;


f = @(x) Powel(x);
gf = @(x) PowelGrad(x);
Hf = @(x) PowelHessian(x);


fbar = 0.001;
c1 = 0.0001;
c2 = 0.1;
stop_tol = 1e-3;

[xopt, fopt, iter] = Newton_95109564(f, gf, Hf, [1 2 2 2], stop_tol, fbar, c1 ,c2)

gradcount
funccount
hesscount


%% SD - Rosen Function
clc
global gradcount;
gradcount = 0;
global funccount;
funccount = 0;
global hesscount;
hesscount = 0;

f = @(x) Rosen(x);
gf = @(x) RoseGrad(x)';


fbar = 0.001;
c1 = 0.0001;
c2 = 0.1;
stop_tol = 1e-3;

[xopt, fopt, iter] = SD_95109564(f, gf, [1,2], stop_tol, fbar, c1 ,c2)

gradcount
funccount
hesscount


%% Newton - Rosen Function
clc
global gradcount;
gradcount = 0;
global funccount;
funccount = 0;
global hesscount;
hesscount = 0;

f = @(x) Rosen(x);
gf = @(x) RoseGrad(x)';
Hf = @(x) RoseHessian(x);


fbar = 0.001;
c1 = 0.0001;
c2 = 0.1;
stop_tol = 1e-3;

[xopt, fopt, iter] = Newton_95109564(f, gf, Hf, [1 2], stop_tol, fbar, c1 ,c2)

gradcount
funccount
hesscount
