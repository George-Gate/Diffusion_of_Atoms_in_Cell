% parameter set 1
clear pars

pars.w=0.2;   % laser beam
pars.D_ph=0;  % no diffusion
pars.T0=20;
pars.P0=0;

initPars;

%% analytical solution
% t,x,y must be scalar
P=@(x,y)pars.T0*pars.P0*exp(-(x.^2+y.^2)*(pars.L/2/pars.w)^2);
anaSol=@(rho_0,t,x,y)expm(-(P(x,y)*pars.matC+pars.alpha*pars.T0*pars.matD)*t)*rho_0 - rho_0;