% parameter set 2 (looks good)
% uniform laser, no diffusion, second type boundary condition with drho_b=0

problemPars=PumpingDiffusionFEMSolver.ProblemPars();
dimRho=problemPars.dimRho;
problemPars.w=200000000;   % uniform laser
problemPars.D_ph=0;  % no diffusion
problemPars.T0=5;
problemPars.boundaryType='second';
problemPars.drho_b_ph=zeros(dimRho,1);

% analytical solution(for no diffusion)
% t,x,y must be scalar
P=@(obj,x,y)obj.T0*obj.P0*exp(-(x^2+y^2)*(obj.L/2/obj.w)^2);
problemPars.anaSolForm=@(obj,t,x,y,z)expm(-(P(obj,x,y)*obj.matC+obj.alpha*obj.T0*obj.matD)*t)*obj.rho_0 - obj.rho_0;


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
