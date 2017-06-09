disp(pwd);
%% physical parameters
parName='parSet2';
%run(['codeVerification/',parName]);
run(parName);

%% simulation parameters
pars.K=10;
% meshGrid={[-1;-pars.w;pars.w;1];[-1;-pars.w;pars.w;1];[-1;1]};  % 9 domains
meshGrid={[-1;1];[-1;1];[-1;1]};  % 1 domain
% meshGrid={[-1;1];[-1;1];[-1;0;1]};  % 5 domains

dimRho=pars.dimRho;
%% initial state and boundary condition
pars.rho_0=ones(dimRho,1)/8/dimRho;

boundaryType='secondZero';
% boundaryType='second';
if strcmp(boundaryType,'firstUniform')
    pars.rho_b=ones(dimRho,1)/8/dimRho;
elseif strcmp(boundaryType,'secondZero')
    % nothing to do
end


%% initialize
% mesh
fprintf('makeMesh...... ');tic;
mesh0=makeMesh3D_cubic(meshGrid{1},meshGrid{2},meshGrid{3});
disp([num2str(toc),' s']);

% coefficient matrix
fprintf('getCoeffs...... ');tic;
% if strcmp(boundaryType,'firstUniform') || strcmp(boundaryType,'first')
    [ M, MM, SS, MG, Nbasis, fun2No, No2fun, getNoByIxyz ]=getCoeffs3D('Bilinear+Lobatto',mesh0,pars.K,dimRho,pars.funP,pars.matC,pars.alpha*pars.T0*pars.matD,boundaryType);
    if ~issparse(MG)
        % use in situ update to save memory
        % H will be full matrix.
        H=MG;
        clear MG
        id=find(SS);
        H(id)=H(id)+pars.D*SS(id);
        H=-H;
        clear id
    else
        % H will be sparse matrix
        % for dense matrix, full matrix is much faster than sparse matrix
        H=-pars.D*SS-MG;
    end
% end
disp([num2str(toc),' s']);

% H=full(H);

% boundary
fprintf('Set boundary...... ');tic;
if strcmp(boundaryType,'firstUniform')
    if sum(abs(pars.rho_0-pars.rho_b))<eps
        ngp=pars.K+10;
        u_0=zeros(Nbasis*dimRho,1);
        Gf=calcGf( mesh0,pars.funP,pars.rho_b,pars.matC,pars.alpha*pars.T0*pars.matD,getNoByIxyz,dimRho,pars.K,Nbasis,ngp );
    else
        error('Unsupported initial condition.');
    end
elseif strcmp(boundaryType,'secondZero')
    Gf=zeros(Nbasis*dimRho,1);
    % expand the initial state by solving linear system M*u=vec{(v_i,rho_0)}
    u_0=zeros(Nbasis*dimRho,1);
    vec1=M\basisInnerProduct( 1, getNoByIxyz, pars.K, Nbasis, mesh0, 3);
    for i=1:dimRho
        u_0(1+(i-1)*Nbasis:i*Nbasis)=pars.rho_0(i)*vec1;
    end
end
disp([num2str(toc),' s']);

%% time evolution
fprintf('Time evolution...... ');tic;
opts=odeset('Mass',MM,'MaxStep',0.1/pars.T0,'InitialStep',1e-6/pars.T0);
ptsPerSecond=50;
t_span=linspace(0,1,ptsPerSecond*pars.T0+1);
[soltList,u] = ode45(@(t,u)H*u-Gf,t_span,u_0,opts);
disp([num2str(toc),' s']);

% shift solution
if strcmp(boundaryType,'secondZero')
    u=u-repmat(u_0.',size(u,1),1);
end

fprintf('Save result to disk...... ');tic;
filename=[parName,' K=',num2str(pars.K),', T0=',num2str(pars.T0),'  boundary=',boundaryType,'.mat'];
save(filename,'pars','meshGrid',...
              'mesh0','boundaryType','Nbasis','fun2No','No2fun','getNoByIxyz',...
              'soltList','u','ptsPerSecond','-v7.3');
disp([num2str(toc),' s']);


%% mapping to space-time domain
dimRho=pars.dimRho;
fprintf('Combine solution...... ');tic;
viewDim='1D';
if strcmp(viewDim,'1D')
    samplingLines={
    %               linspace(-1,1,200),zeros(1,200),zeros(1,200),'xAxis';
                  zeros(1,200),linspace(-1,1,200),zeros(1,200),'yAxis';
    %               zeros(1,200),zeros(1,200),linspace(-1,1,200),'zAxis';
    %               linspace(-1,1,200),zeros(1,200),-0.99*ones(1,200),'xAxis2';
%                    0,0,0,'centerPoint';
%                    0.5,0.5,0.5,'(0.5,0.5,0.5)';
%                    1,1,1,'(1,1,1)';
     };
    stmax=pars.T0;
    sNpoints=100;
    tid=[1:max(1,round(ptsPerSecond*stmax/sNpoints)):ptsPerSecond*stmax,ptsPerSecond*stmax+1];
    sol=combineSolution( u(tid,:), '1D', mesh0, getNoByIxyz, pars.K, Nbasis, dimRho, samplingLines);
elseif strcmp(viewDim,'3D')
    tid=125;
    sol=combineSolution( u(tid,:), '3D', mesh0, getNoByIxyz, pars.K, Nbasis, dimRho, 50 );
end
combinedTime=reshape(soltList(tid),length(tid),1);
disp([num2str(toc),' s']);

%% plot
resultPlotter( sol, pars, combinedTime );

%% get analytical solution
if strcmp(sol.sampleDim,'1D')
    anaSol_diff.sampleDim='1D';
    for Lid=1:length(sol.tag)
        px=sol.xList{Lid};py=sol.yList{Lid};pz=sol.zList{Lid};
        anaSol_diff.tag{Lid}=[sol.tag{Lid},'(diff with anaSol)'];
        anaSol_diff.xList{Lid}=px;
        anaSol_diff.yList{Lid}=py;
        anaSol_diff.zList{Lid}=pz;
        anaSol_diff.rho{Lid}=zeros(length(px),length(combinedTime),pars.dimRho);
        for kkk=1:length(combinedTime)
            for ipx=1:length(px)
                anaSol_diff.rho{Lid}(ipx,kkk,:)=reshape(anaSol(pars.rho_0,combinedTime(kkk),px(ipx),py(ipx)),1,1,pars.dimRho)-sol.rho{Lid}(ipx,kkk,:);
    %            anaSol_diff.rho{Lid}(ipx,kkk,:)=reshape(anaSol(pars.rho_0,combinedTime(kkk),px(ipx),py(ipx)),1,1,pars.dimRho);
            end
        end
    end
    resultPlotter( anaSol_diff, pars, combinedTime );
    close;
end