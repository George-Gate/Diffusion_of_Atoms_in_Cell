% parameter set 1 
% second type boundary condition with drho_b=0

problemPars=PumpingDiffusionFEMSolver.ProblemPars();
dimRho=problemPars.dimRho;
problemPars.w=0.2;   % uniform laser
problemPars.D_ph=0.1;  % no diffusion
problemPars.T0=5;
problemPars.L=1;
problemPars.boundaryType='second';
problemPars.rho_b_ph=zeros(dimRho,1);
problemPars.rho_0_ph=ones(dimRho,1)/dimRho/problemPars.L^3;
