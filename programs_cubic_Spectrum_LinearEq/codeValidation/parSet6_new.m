% parameter set 6
% gaussian laser, diffusion, robin type boundary condition with 
%   drho/dn=-A*( rho-ones(dimRho,1) )

problemPars=PumpingDiffusionFEMSolver.ProblemPars();
dimRho=problemPars.dimRho;
problemPars.w=0.8;   % gaussian laser
problemPars.D_ph=0.1;  % diffusion
problemPars.T0=5;
problemPars.boundaryType='robin';
problemPars.robinF_ph=ones(dimRho,1)/dimRho;
problemPars.robinA_ph=0.5;



%% plot anaSol
% dim=2;
% [xList,yList]=meshgrid(-1:0.01:1);
% anaSol_num=zeros(size(xList));
% for i=1:length(xList)
%     for j=1:length(yList)
%         tmp=problemPars.anaSol(2,xList(j,i),yList(j,i),0);
%         anaSol_num(j,i)=tmp(dim);
%     end
% end
% %%
% figure();
% mesh(xList,yList,anaSol_num);
