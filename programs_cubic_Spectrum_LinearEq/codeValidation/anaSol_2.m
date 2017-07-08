% parameter set for analytical solution 2
% see .docx file for details


D_ph=0.1;
problemPars=PumpingDiffusionFEMSolver.ProblemPars();
% set G matrix
problemPars.matC=zeros(2);
problemPars.matD=[4,-14*D_ph;-4,14*D_ph];
problemPars.alpha=1;
problemPars.P0=0;

problemPars.D_ph=D_ph;  % diffusion
problemPars.T0=5;
problemPars.L=1.8;
problemPars.boundaryType='robin';
problemPars.robinA_ph={2;
                       [0,2;2,0];
                       3;
                       -3;
                       1;
                       -1};
problemPars.robinF_ph=zeros(2,1);
problemPars.rho_0_ph={{@(x)-exp(2*x),@(y)exp(3*y),@(z)exp(z)};
                      {@(x)exp(2*x),@(y)exp(3*y),@(z)exp(z)}};

% analytical solution
problemPars.anaSolForm=@(obj,t,x,y,z)[-exp(2*x+3*y+z-4*t);     % x,y,z and t should be row vectors
                                       exp(2*x+3*y+z-4*t)];


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
