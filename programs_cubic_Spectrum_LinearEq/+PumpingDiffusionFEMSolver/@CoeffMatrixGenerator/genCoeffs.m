function [  ] = genCoeffs( obj )
%genCoeffs  To generate all coefficient matrices
    
    % generate mesh
    if (obj.mesh.outdated)
        obj.mesh.makeMesh();
    end
    % base numbering
    if ~(obj.simuPars.maxOrder==obj.baseFunHandle.maxOrder ...
        && obj.mesh.meshPars_current==obj.baseFunHandle.meshPars)
        Nbasis=obj.baseFunHandle.numbering(obj.mesh,obj.simuPars.maxOrder);
    end

    % ----------------------- Gen MM, SS, MG ---------------------------------------
    startT=tic;
    % call genCoeffs_boundary_independent to get M, S, MP
    % use mex-function? No need!
    obj.genCoeffs_boundary_independent();
    
    % matrix statistics
    import PumpingDiffusionFEMSolver.library.matStat
    matStat(obj.M,'M');matStat(obj.S,'S');matStat(obj.MP,'MP');
    
    % calc MM and SS
    obj.MM=obj.M; obj.SS=obj.S;
    for i=1:obj.problemPars.dimRho-1
        obj.MM=blkdiag(obj.MM,obj.M);
        obj.SS=blkdiag(obj.SS,obj.S);
    end
    
    % display time consumption
    import PumpingDiffusionFEMSolver.library.sec2hms;
    disp(['Time used to generate matrix MM, SS and MP: ',sec2hms(toc(startT))]);
    
    startT=tic;
    % calc MG 
    matDensity=nnz(obj.MP)/numel(obj.MP);
    D0=obj.problemPars.alpha * obj.problemPars.T0 * obj.problemPars.matD;
    C0=obj.problemPars.matC;
    if (matDensity*0.5938>0.2)
        % use full matrix to save memory
        disp('Use full matrix for MG.');
        MG=kron(full(C0),full(obj.MP));
        MG_2=kron(sparse(D0),obj.M);
        id=find(MG_2);
        MG(id)=MG(id)+MG_2(id);
        clear id MG_2
    else
        % kron() for sparse matrix needs x5 extra memory for intermediate variables
        % x1 means the size of non-zero elements, the size of sparse matrix is x2.
        disp('Use sparse matrix for MG.');
        MG=kron(sparse(C0),sparse(obj.MP));
        MG=MG+kron(sparse(D0),obj.M);
    end
    obj.MG=MG;
    clear MG
    
    % matrix statistics
    matStat(obj.MG,'MG');
    
    % display time consumption
    disp(['Time used to generate matrix MG: ',sec2hms(toc(startT))]);
    
    % ---------------------- Gen boundary term ------------------------------------
    % call genVecQ/genCBvecR/genMABvecF according to boundaryType
    switch obj.problemPars.boundaryType
        case 'first'
            startT=tic;
            obj.genCBvecR();
            disp(['Time used to generate CB and vecR: ',sec2hms(toc(startT))]);
            matStat(obj.vecR,'vecR');matStat(obj.CB,'CB');
        case 'second'
            startT=tic;
            obj.genVecQ();
            disp(['Time used to generate vecQ: ',sec2hms(toc(startT))]);
            matStat(obj.vecQ,'vecQ');
        case 'robin'
            startT=tic;
            obj.genMABvecF();
            disp(['Time used to generate MAB and vecF: ',sec2hms(toc(startT))]);
            matStat(obj.MAB,'MAB');matStat(obj.vecF,'vecF');
        otherwise
            error(['Unknow boundary type: ',obj.problemPars.boundaryType]);
    end
    

    % ---------------------- Record current pars -----------------------------------
    obj.simuPars_current=obj.simuPars;
    obj.problemPars_current=obj.problemPars;
    obj.meshPars_current=obj.mesh.meshPars_current;
end

