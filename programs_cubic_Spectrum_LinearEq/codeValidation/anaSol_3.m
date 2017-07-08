% parameter set for analytical solution 2
% see .docx file for details


D_ph=0.1;
L=1.8;
problemPars=PumpingDiffusionFEMSolver.ProblemPars();
% set G matrix
problemPars.matC=[0,0,-1;0,0,1;0,0,0];
problemPars.matD=[1,0,0;-1,1,0;0,0,0];
problemPars.funPForm_ph={@(obj,x)1,@(obj,y)1,@(obj,z)1./z};
problemPars.alpha=1;
problemPars.P0=1;  % enable matrix MP calc


problemPars.D_ph=D_ph;  % diffusion
problemPars.T0=5;
problemPars.L=L;
problemPars.boundaryType='second';
problemPars.drho_b_ph={{ {@(y)y,@(z)-1} ; 0 ; {@(y)y,@(z)-z} };   % [-y;0;-yz]
                       { {@(y)y,@(z) 1} ; 0 ; {@(y)y,@(z) z} };   % [ y;0; yz]
                       { {@(z)-1,@(x)x} ; 0 ; {@(z)-z,@(x)x} };   % [-x;0;-xz]
                       { {@(z) 1,@(x)x} ; 0 ; {@(z) z,@(x)x} };   % [ x;0; xz]
                       { 0 ; 0 ; {@(x)-x,@(y)y} };                % [ 0;0;-xy]
                       { 0 ; 0 ; {@(x) x,@(y)y} } };              % [ 0;0; xy]

rho0_1={ {@(x)x,@(y)y,@(z)1};  0; {@(x)x,@(y)y,@(z)z} };
rho0_2=[1;0;0];
problemPars.rho_0_ph={rho0_1,rho0_2};   % rho_0_ph=rho0_1+rho0_2=[xy;0;xyz]+[1;0;0]

% analytical solution
problemPars.anaSolForm=@(obj,t,x,y,z)[x.*y+exp(-t);     % x,y,z and t should be row vectors
                                      t.*exp(-t)+zeros(size(x));
                                      x.*y.*z+zeros(size(t))];


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
