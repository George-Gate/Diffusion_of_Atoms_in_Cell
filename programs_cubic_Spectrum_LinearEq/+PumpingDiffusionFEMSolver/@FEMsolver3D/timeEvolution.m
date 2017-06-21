function [ ] = timeEvolution( obj,filename, u0 )
%timeEvolution  Perform the time evolution
    import PumpingDiffusionFEMSolver.library.sec2hms;
    
    % makeMesh, baseNumbering, generateCoeffMatrices
    % ##!! POSSIBLE BUG: If someone change the baseNumbering but not other parameters, then
    %      coeffMatrix will not be outdated, thus may cause errors. DO NOT use the numbering()
    %      method of @BaseFunction class manually.  #<FIXED>
    if obj.coeffMatrix.outdated
        obj.coeffMatrix.genCoeffs();
    elseif ~(obj.coeffMatrix.simuPars_current.maxOrder==obj.baseFunction.maxOrder ...     % in case that someone change the baseNumbering but not other parameters
        && obj.coeffMatrix.meshPars_current==obj.baseFunction.meshPars)
        obj.coeffMatrix.genCoeffs();
    end
    
    % get initial vector
    if exist('u0','var') && isnumeric(u0) && iscolumn(u0) ...
             && length(u0)==size(obj.MM,1)
        obj.u0=u0;
    else
        startT=tic;
        obj.u0=obj.getInitialState();
        disp(['Time used to calc initial state: ',sec2hms(toc(startT))]);
    end
    
    
    % EQ: du/dt=Hu+F
    startT=tic;
    pars=obj.problemPars;
    switch obj.problemPars.boundaryType
        case 'first'
            obj.H=-pars.D*obj.SS-obj.MG+pars.D*obj.CB;
            obj.F=-obj.vecR;
        case 'second'
            obj.H=-pars.D*obj.SS-obj.MG;
            obj.F=obj.vecQ;
        case 'robin'
            obj.H=-pars.D*obj.SS-obj.MG-pars.D*obj.MAB;
            obj.F=obj.vecF;
        otherwise
            error(['Unknow boundary type ',obj.problemPars.boundaryType]);
    end
    disp(['Time used to set ODE systems: ',sec2hms(toc(startT))]);

    % time evolution (H and F should be time independent)
    startT=tic;
    disp('Time evolution started.');
    opts=odeset('Mass',obj.MM,'MaxStep',0.1/pars.T0,'InitialStep',1e-6/pars.T0);
    t_span=linspace(0,1,obj.sampleRate*pars.T0+1);
    [obj.sol_t,obj.sol_u] = ode45(@(t,u)obj.H*u+obj.F,t_span,obj.u0,opts);
    disp(['Time used for evolution: ',sec2hms(toc(startT))]);
    
    % save result to file
    startT=tic;
    if ~exist('filename','var') || ~isa(filename,'char') || strcmp(filename,'')
        filename=['FEMResult_',datestr(datetime,'yyyymmdd_HH_MM_SS'),'.mat'];
    end
    FEMResult.sol_t=obj.sol_t;
    FEMResult.sol_u=obj.sol_u;
    FEMResult.u0=obj.u0;
    FEMResult.meshPars=obj.meshPars;
    FEMResult.simuPars=obj.simuPars;
    FEMResult.simuPars.baseFunHandle=obj.simuPars.baseFunHandle.copy;  % use deep copy to avoid cross affecting
    FEMResult.problemPars=obj.problemPars;
    FEMResult.getNoByIxyz=obj.baseFunction.getNoByIxyz;
    FEMResult.Nbasis=obj.baseFunction.Nbasis;
    FEMResult.sampleRate=obj.sampleRate;
    save(filename,'FEMResult');
    obj.FEMResult=FEMResult;
    disp(['Evolution result and parameters saved to "',filename,'"']);
    disp(['Time used to save result: ',sec2hms(toc(startT))]);
    
end

