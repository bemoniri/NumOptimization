clc
clear all
format long
addpath('./functions')

%% BFGS - Powel Function

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
stop_tol = 0.1;

[xopt, fopt, iter] = BFGS_95109564(f, gf, [1;2;2;2], stop_tol,  fbar, c1 ,c2)

gradcount
funccount
hesscount

%% BFGS - Rosen Function

global gradcount;
gradcount = 0;
global funccount;
funccount = 0;
global hesscount;
hesscount = 0;

f = @(x) Rosen(x);
gf = @(x) RoseGrad(x);

fbar = 0.01;
c1 = 0.0001;
c2 = 0.1;
stop_tol = 0.1;

[xopt, fopt, iter] = BFGS_95109564(f, gf, [1;2], stop_tol,  fbar, c1 ,c2)

gradcount
funccount
hesscount
