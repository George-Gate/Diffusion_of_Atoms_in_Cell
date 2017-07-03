% parameter set 5
% gaussian laser, diffusion, robin type boundary condition with 
%   drho_i/dn=-sum_k( Aik*(rho_i-rho_k) )

problemPars=PumpingDiffusionFEMSolver.ProblemPars();
dimRho=problemPars.dimRho;
problemPars.w=0.2;   % gaussian laser
problemPars.D_ph=0.1;  % diffusion
problemPars.T0=5;
problemPars.boundaryType='robin';
problemPars.robinF=zeros(dimRho,1);
A=zeros(dimRho,dimRho);
channel=[5,2,0.1;   % 5 <-> 2, 0.1
         5,8,0.3];  % 5 <-> 8, 0.3
for iii=1:size(channel,1)
    i=channel(iii,1);k=channel(iii,2);rate=channel(iii,3);
    A(i,k)=A(i,k)-rate;
    A(i,i)=A(i,i)+rate;
    A(k,i)=A(k,i)-rate;
    A(k,k)=A(k,k)+rate;
end
problemPars.robinA=A;



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
