classdef CoeffMatrixGenerator < handle
    %CoeffMatrixGenerator  Generate coefficient matrices for FEMsolver3D
    %   
    
    properties(SetAccess=private)
        MM=[]; SS=[]; MG=[];   % coeff matrix
        vecQ=[];           % For second type boundary condition
        CB=[];  vecR=[];   % For first type boundary condition
        MAB=[]; vecF=[];   % For Robin type boundary condition
        M=[]; S=[]; MP=[]; % Intermediate variables
        simuPars_current=0;    % the simuPars for currently generated coeff matrices
        problemPars_current=0; % the problemPars for currently generated coeff matrices
        meshPars_current=0;    % the meshPars for currently generated coeff matrices
    end
    
    properties
        simuPars;       % parameters for simulation, class @SimuPars
        problemPars;    % parameters for the problem, class @ProblemPars
        mesh;           % @Mesh3D object.
    end
    
    properties(Dependent)
        outdated;   % return whether the currently generated coeff matrices are outdated.
        baseFunHandle;  % return the function hande in simuPars
    end
    
    methods
        function obj=CoeffMatrixGenerator(simuPars, problemPars, mesh)
            if nargin==3
                obj.simuPars=simuPars;
                obj.problemPars=problemPars;
                obj.mesh=mesh;
            else
                obj.simuPars=PumpingDiffusionFEMSolver.SimuPars();
                obj.problemPars=PumpingDiffusionFEMSolver.ProblemPars();
                obj.mesh=PumpingDiffusionFEMSolver.Mesh3D();
            end
        end
        function set.simuPars(obj,val)
            if isa(val,'PumpingDiffusionFEMSolver.SimuPars')
                obj.simuPars=val;
            else
                error('simuPars should be an object of class @PumpingDiffusionFEMSolver.SimuPars');
            end
        end
        function set.problemPars(obj,val)
            if isa(val,'PumpingDiffusionFEMSolver.ProblemPars')
                obj.problemPars=val;
            else
                error('simuPars should be an object of class @PumpingDiffusionFEMSolver.ProblemPars');
            end
        end
        function set.mesh(obj,val)
            if isa(val,'PumpingDiffusionFEMSolver.Mesh3D')
                obj.mesh=val;
            else
                error('mesh should be an object of class @PumpingDiffusionFEMSolver.Mesh3D');
            end
        end
        function outdated=get.outdated(obj)
            outdated=~( obj.simuPars==obj.simuPars_current && obj.problemPars==obj.problemPars_current && obj.mesh.meshPars_current==obj.meshPars_current);
        end
        function fun=get.baseFunHandle(obj)
            fun=obj.simuPars.baseFunHandle;
        end
        
        % To generate all coefficient matrices
        genCoeffs(obj);
    end
    
    methods(Hidden, Access=private)
        genCoeffs_boundary_independent(obj);
        genCoeffs_boundary_independent_slow(obj);
        % The following methods will be called by genCoeffs() method automatically according to the boundary type.
        genVecQ(obj);     % generate coefficient matrices for second type boundary condition
        genCBvecR(obj);   % generate coefficient matrices for first type boundary condition
        genMABvecF(obj);  % generate coefficient matrices for Robin type boundary condition
    end
    
end

