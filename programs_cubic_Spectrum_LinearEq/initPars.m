% set default parameters and load matC and matD

load DiffusionFEM.mat;  % load matC and matD

%% default physical parameters
pars.matC=matC;
pars.matD=matD;
if ~isfield(pars,'alpha')
    pars.alpha=0.5; % 0.5/s
end
if ~isfield(pars,'w')
    pars.w=0.2;   % 0.2cm
end
if ~isfield(pars,'P0')
    pars.P0=1;    % 1/s
end
if ~isfield(pars,'D_ph')
    pars.D_ph=0.1; % 0.1 cm^2/s
end

pars.dimRho=size(matD,1);

if ~isfield(pars,'T0')
    pars.T0=0.5;   % evolution time, in unit of seconds
end
if ~isfield(pars,'L')
    pars.L=1;     % length of cubic cell, in unit of cm
end

pars.funP={@(x)exp( -x.*x * (pars.L/2/pars.w)^2 );...
           @(y)exp( -y.*y * (pars.L/2/pars.w)^2 );...
           @(z)pars.T0*pars.P0};

% update pars.D
pars.D=4*pars.T0/pars.L/pars.L*pars.D_ph;

% pars.alpha*pars.T0*matD -> getCoeffs3D

clear matC matD;