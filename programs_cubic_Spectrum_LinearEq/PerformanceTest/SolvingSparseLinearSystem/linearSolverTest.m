parSet1;

FEMSolver=PumpingDiffusionFEMSolver.FEMsolver3D(problemPars);
FEMSolver.mesh.use_mex=0;
FEMSolver.mesh.meshPars.Nx=2;
FEMSolver.mesh.meshPars.Ny=1;
FEMSolver.mesh.meshPars.Nz=1;
FEMSolver.simuPars.maxOrder=16;
FEMSolver.genCoeffs();
u0=FEMSolver.getInitialState();
pars=FEMSolver.problemPars;
switch FEMSolver.problemPars.boundaryType
    case 'first'
        H=-pars.D*FEMSolver.SS-FEMSolver.MG+pars.D*FEMSolver.CB;
        F=-FEMSolver.vecR;
    case 'second'
        H=-pars.D*FEMSolver.SS-FEMSolver.MG;
        F=FEMSolver.vecQ;
    case 'robin'
        H=-pars.D*FEMSolver.SS-FEMSolver.MG-pars.D*FEMSolver.MAB;
        F=FEMSolver.vecF;
    otherwise
        error(['Unknow boundary type ',FEMSolver.problemPars.boundaryType]);
end
MM=FEMSolver.MM;
M=FEMSolver.M;
disp('Coeff calc finished.');

%%
tol=1E-6;

%% direct solver
tID=1;
name{tID}='Direct Solver';
t{tID}=timeit(@()MM\(H*u0-F),1);
u{tID}=MM\(H*u0-F);
err{tID}=max(abs( MM*u{tID} - (H*u0-F) ))/norm(u{tID});
disp(['t',num2str(tID),'=',num2str(t{tID}),'s  err',num2str(tID),'=',num2str(err{tID})]);
iter{tID}=0;

%% direct solver + LU
tID=2;
name{tID}='Direct Solver + LU';
[L,U,P,Q,R] = lu(MM);
t{tID}=timeit(@()Q * ( U \ (L \ (P * (R \(H*u0-F) )))),1);
u{tID}=Q * ( U \ (L \ (P * (R \(H*u0-F) ))));
err{tID}=max(abs( MM*u{tID} - (H*u0-F) ))/norm(u{tID});
disp(['t',num2str(tID),'=',num2str(t{tID}),'s  err',num2str(tID),'=',num2str(err{tID})]);
iter{tID}=0;

%% pcg without preconditioner
tID=3;
name{tID}='pcg';
u00=u0+u{1}*0.001;
t{tID}=timeit(@()pcg(MM,(H*u00-F),tol,1000,[],[],u{1}),2);
[u{tID},flag{tID},relres{tID},iter{tID},resvec{tID}]=pcg(MM,(H*u00-F),tol,1000,[],[],u{1});
err{tID}=max(abs( MM*u{tID} - (H*u00-F) ))/norm(u{tID});
disp(['t',num2str(tID),'=',num2str(t{tID}),'s  err',num2str(tID),'=',num2str(err{tID})]);

%% pcg + incomplete Cholesky factorization with zero-fill as preconditioner
tID=4;
name{tID}='pcg + ichol';
u00=u0+u{1}*0.001;
L = ichol(MM);
t{tID}=timeit(@()pcg(MM,(H*u00-F),tol,1000,L,L',u{1}),2);
[u{tID},flag{tID},relres{tID},iter{tID},resvec{tID}]=pcg(MM,(H*u00-F),tol,1000,L,L',u{1});
err{tID}=max(abs( MM*u{tID} - (H*u00-F) ))/norm(u{tID});
disp(['t',num2str(tID),'=',num2str(t{tID}),'s  err',num2str(tID),'=',num2str(err{tID})]);


%% pcg + ict incomplete Cholesky preconditioner as preconditioner
tID=5;
name{tID}='pcg + ichol(ict)';
u00=u0+u{1}*0.001;
L2 = ichol(MM,struct('type','ict','droptol',0.001,'diagcomp',2));
t{tID}=timeit(@()pcg(MM,(H*u00-F),tol,1000,L2,L2',u{1}),2);
[u{tID},flag{tID},relres{tID},iter{tID},resvec{tID}]=pcg(MM,(H*u00-F),tol,1000,L2,L2',u{1});
err{tID}=max(abs( MM*u{tID} - (H*u00-F) ))/norm(u{tID});
disp(['t',num2str(tID),'=',num2str(t{tID}),'s  err',num2str(tID),'=',num2str(err{tID})]);

%% minres without preconditioner
tID=6;
name{tID}='minres';
u00=u0+u{1}*0.001;
t{tID}=timeit(@()minres(MM,(H*u00-F),tol,1000,[],[],u{1}),2);
[u{tID},flag{tID},relres{tID},iter{tID},resvec{tID}]=minres(MM,(H*u00-F),tol,1000,[],[],u{1});
err{tID}=max(abs( MM*u{tID} - (H*u00-F) ))/norm(u{tID});
disp(['t',num2str(tID),'=',num2str(t{tID}),'s  err',num2str(tID),'=',num2str(err{tID})]);


%% minres + incomplete Cholesky factorization with zero-fill as preconditioner
tID=7;
name{tID}='minres + ichol';
u00=u0+u{1}*0.001;
L3 = ichol(MM);
t{tID}=timeit(@()minres(MM,(H*u00-F),tol,1000,L3,L3',u{1}),2);
[u{tID},flag{tID},relres{tID},iter{tID},resvec{tID}]=pcg(MM,(H*u00-F),tol,1000,L3,L3',u{1});
err{tID}=max(abs( MM*u{tID} - (H*u00-F) ))/norm(u{tID});
disp(['t',num2str(tID),'=',num2str(t{tID}),'s  err',num2str(tID),'=',num2str(err{tID})]);


%% symmlq without preconditioner
tID=8;
name{tID}='symmlq';
u00=u0+u{1}*0.001;
t{tID}=timeit(@()symmlq(MM,(H*u00-F),tol,1000,[],[],u{1}),2);
[u{tID},flag{tID},relres{tID},iter{tID},resvec{tID}]=symmlq(MM,(H*u00-F),tol,1000,[],[],u{1});
err{tID}=max(abs( MM*u{tID} - (H*u00-F) ))/norm(u{tID});
disp(['t',num2str(tID),'=',num2str(t{tID}),'s  err',num2str(tID),'=',num2str(err{tID})]);

%% symmlq with diag preconditioner
tID=9;
name{tID}='symmlq + diag(diag(MM))';
u00=u0+u{1}*0.001;
L4 = (diag(diag(MM)));
t{tID}=timeit(@()symmlq(MM,(H*u00-F),tol,1000,L4,[],u{1}),2);
[u{tID},flag{tID},relres{tID},iter{tID},resvec{tID}]=symmlq(MM,(H*u00-F),tol,1000,L4,[],u{1});
err{tID}=max(abs( MM*u{tID} - (H*u00-F) ))/norm(u{tID});
disp(['t',num2str(tID),'=',num2str(t{tID}),'s  err',num2str(tID),'=',num2str(err{tID})]);

%% output summary
Nx=FEMSolver.mesh.meshPars.Nx;Ny=FEMSolver.mesh.meshPars.Ny;Nz=FEMSolver.mesh.meshPars.Nz;
K=FEMSolver.simuPars.maxOrder;
disp('========================[Test Summary]===========================');
!hostname > hostname
hostname=textread('hostname','%s'); %#ok<DTXTRD>
delete hostname
fprintf('Test host: %s\n',hostname{1});
fprintf('Basis Configuration: %dx%dx%d, maxOrder=%d\n',Nx,Ny,Nz,K);
fprintf('tol=%G\n',tol);
matStat(MM);matStat(H);
fprintf('-----------------------------------------------------\n')
fprintf(' ID |  Time  |    Err   | N iter. |   Description              \n');
fprintf('----+--------+----------+---------+------------------\n');
for tID=1:length(t)
fprintf('%3d | %6.2f | %8.2E | %7d | %s\n',tID,t{tID},err{tID},iter{tID},name{tID});
end
disp('=======================[End of Summary]==========================');