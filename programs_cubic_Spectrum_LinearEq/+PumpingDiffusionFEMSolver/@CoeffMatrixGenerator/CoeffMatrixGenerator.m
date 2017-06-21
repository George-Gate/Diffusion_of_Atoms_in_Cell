classdef CoeffMatrixGenerator < handle
    %CoeffMatrixGenerator  Generate coefficient matrices for FEMsolver3D
    
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
        % ###TODO: need to define the eq operator for @SimuPars and @ProblemPars
        outdated;   % return whether the currently generated coeff matrices are outdated.
        baseFunHandle;  % return the function hande in simuPars
        dimRho;
    end
    
    methods
        function obj=CoeffMatrixGenerator(varargin)
            if nargin>3
                error('Too many input arguments, maximum is three.');
            end
            meshSet=false;simuSet=false;probSet=false;
            for i=1:nargin
                if isa(varargin{i},'PumpingDiffusionFEMSolver.MeshPars') && ~meshSet
                    obj.mesh=PumpingDiffusionFEMSolver.Mesh3D(varargin{i});
                    meshSet=true;
                elseif isa(varargin{i},'PumpingDiffusionFEMSolver.Mesh3D') && ~meshSet
                    obj.mesh=varargin{i};
                    meshSet=true;
                elseif isa(varargin{i},'PumpingDiffusionFEMSolver.SimuPars') && ~simuSet
                    obj.simuPars=varargin{i};
                    simuSet=true;
                elseif isa(varargin{i},'PumpingDiffusionFEMSolver.ProblemPars') && ~probSet
                    obj.problemPars=varargin{i};
                    probSet=true;
                else
                    error(['Inputs should be objects of the following class: @MeshPars, @Mesh3D, @SimuPars, @ProblemPars. ', ...
                           'MeshPars and Mesh3D can not be specified at the same time.']);
                end
            end
            if ~meshSet
                obj.mesh=PumpingDiffusionFEMSolver.Mesh3D();
            end
            if ~simuSet
                obj.simuPars=PumpingDiffusionFEMSolver.SimuPars();
            end
            if ~probSet
                obj.problemPars=PumpingDiffusionFEMSolver.ProblemPars();
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
            outdated=~( obj.simuPars==obj.simuPars_current ...
                     && obj.problemPars==obj.problemPars_current ...
                     && obj.mesh.meshPars_current==obj.meshPars_current ...
                     && ~obj.mesh.outdated);
        end
        function fun=get.baseFunHandle(obj)
            fun=obj.simuPars.baseFunHandle;
        end
        function dimRho=get.dimRho(obj)
            dimRho=obj.problemPars.dimRho;
        end
        
        
        % clean coefficient matrices to free memory
        function cleanCoeffMatrices(obj)
            obj.S=[];obj.M=[];obj.SS=[];obj.MM=[];obj.MP=[];obj.MG=[];
            obj.vecQ=[];obj.vecR=[];obj.vecF=[];
            obj.CB=[];obj.MAB=[];
            % set current parameter
            obj.simuPars_current=0;    
            obj.problemPars_current=0; 
            obj.meshPars_current=0;
        end
        
        
        
        % To generate all coefficient matrices
        genCoeffs(obj);
        % To calc the vecX's in boundary terms
        [ vecX ] = baseProjectionOnBoundary( obj, funX, baseFunEvaluator, baseDerivative  )
        
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

