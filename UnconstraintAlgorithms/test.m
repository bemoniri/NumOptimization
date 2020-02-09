clc
clear all
format long

%% SD - Rosen Function
clc
global gradcount;
gradcount = 0;
global funccount;
funccount = 0;
global hesscount;
hesscount = 0;

f = @(x) Rosen(x);
gf = @(x) RoseGrad(x);

[xopt, fopt, iter] = SD(f, gf, [1,2], 1e-3, 1e-5)
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
gf = @(x) RoseGrad(x);
Hf = @(x) RoseHessian(x);

[xopt, fopt, iter] = Newton(f, gf, Hf, [1 2], 10e-3, 10e-5)
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

f = @(x) Powel(x);
gf = @(x) PowelGrad(x)';


[xopt, fopt, iter] = SD(f, gf, [1,2,2,2], 1e-3, 1e-5)
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
gf = @(x) PowelGrad(x)';
Hf = @(x) PowelHessian(x);

[xopt, fopt, iter] = Newton(f, gf,Hf,  [1,2,2,2], 1e-3, 1e-5)
gradcount
funccount
hesscount
