load DiffusionFEM.mat;  % load matC and matD

%% physical parameters
alpha=0.5; % 0.5/s
w=0.02;   % 0.2cm
P0=1;    % 1/s
D_ph=0.1; % 0.1 cm^2/s
dimRho=length(matD);

T0=10;   % evolution time, in unit of seconds
L=1;     % length of cubic cell, in unit of cm

funP={@(x)exp( -x.*x * (L/2/w)^2 );...
      @(y)exp( -y.*y * (L/2/w)^2 );...
      @(z)T0*P0};
D=4*T0/L/L*D_ph;

% alpha*T0*matD -> getCoeffs3D

%% simulation parameters
K=10;
% meshGrid={[-1;-w;w;1];[-1;-w;w;1];[-1;1]};  % 9 domains
meshGrid={[-1;1];[-1;1];[-1;1]};  % 1 domain
% meshGrid={[-1;1];[-1;1];[-1;0;0.5;0.6;0.7;1]};  % 5 domains


%% initial state and boundary condition
boundaryType='firstUniform';
% boundaryType='second';
rho_0=ones(dimRho,1)/8/dimRho;
rho_b=ones(dimRho,1)/8/dimRho;


%% initialize
% mesh
tic;
mesh0=makeMesh3D_cubic(meshGrid{1},meshGrid{2},meshGrid{3});
disp(['Time to makeMesh: ',num2str(toc),' s']);tic;
% coefficient matrix
% if strcmp(boundaryType,'firstUniform') || strcmp(boundaryType,'first')
    [ MM, SS, MG, Nbasis, fun2No, No2fun, getNoByIxyz_a ]=getCoeffs3D('Bilinear+Lobatto',mesh0,K,dimRho,funP,matC,alpha*T0*matD,boundaryType);
    H=-D*SS-MG;
% end
disp(['Time to getCoeffs: ',num2str(toc),' s']);

tic;
% boundary
if strcmp(boundaryType,'firstUniform')
    if sum(abs(rho_0-rho_b))<eps
        ngp=K+10;
        u_0=zeros(Nbasis*dimRho,1);
        Gf=calcGf( mesh0,funP,rho_b,matC,alpha*T0*matD,getNoByIxyz_a,dimRho,K,Nbasis,ngp );
    else
        error('Unsupported initial condition.');
    end
end
disp(['Time to Set Boundary: ',num2str(toc),' s']);
%% time evolution
opts=odeset('Mass',MM);
tic;
[tList,u]=ode45(@(t,u)H*u-Gf,[0,1],u_0,opts);
disp(['Time to calc evolution: ',num2str(toc),' s']);

tic;
filename=['K=',num2str(K),'.mat'];
save(filename,'alpha','w','P0','D_ph','dimRho','T0','L','funP','K','meshGrid',...
              'mesh0','boundaryType','rho_0','rho_b','Nbasis','fun2No','No2fun','getNoByIxyz_a',...
              'tList','u');
disp(['Time to save result to disk: ',num2str(toc),' s']);

%% mapping to space-time domain
tic;
tid=[10000];
combinedTime=reshape(tList(tid),length(tid),1);
sol=combineSolution( u(tid,:), '3D', mesh0, getNoByIxyz_a, K, Nbasis, dimRho, 100 );
disp(['Time to combine solution: ',num2str(toc),' s']);

%% plot
hf=figure('position',[1426 82 1811 533]);
tid=1;
for rhoComponent=1:8
    subplot(2,4,rhoComponent);

    for Did=1:mesh0.Ndomains
        midPos=sum(mesh0.domains.xyz(:,:,Did),2)/2;
        domW=mesh0.domains.h(:,Did)/2;
        xGrid=sol.xSample*domW(1)+midPos(1);
        yGrid=sol.ySample*domW(2)+midPos(2);
        zGrid=sol.zSample*domW(3)+midPos(3);
        hs=slice(xGrid,yGrid,zGrid,sum(sol.rho{Did}(:,:,:,tid,rhoComponent),5),midPos(1),midPos(2),midPos(3));
        for i=1:length(hs)
            hs(i).EdgeColor='none';
        end
    end
    view([26 6]);
    colormap winter
    colorbar;
    title(['\rho_',num2str(rhoComponent)]);
    xlabel('x');ylabel('y');zlabel('z');
end
titleTextBox=uicontrol(hf,'style','text');
annotation('textbox','String',{['\rho(t)-\rho(0),  K=',num2str(K)],['T=',num2str(T0*combinedTime(tid)),' s']},'Position',[0.02 0.8 0.1 0.1]);