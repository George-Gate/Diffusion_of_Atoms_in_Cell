classdef FEMsolver3D < handle
    %FEMsolver3D   FEM simulator
    
% --------------- Version Check ---------------------------
    properties(Constant)
        classVersion='0.1';
    end
    properties(SetAccess=private)
        objectVersion='0.1';
    end
    
    methods
        function set.objectVersion(obj,val)
            if ~strcmp(val,obj.classVersion)
                warning('objectVersion ~= classVersion, the definition of class may changed.');
            end
            obj.objectVersion=val;
        end
    end
% --------------- End of Version Check -------------------

    properties(Access=private)
        coeffMatrix;    % @CoeffMatrixGenerator object
    end
    
    properties(SetAccess=private)
        H;
        F;
        u0;
        sol_u;
        sol_t;
        FEMResult;
    end
    
    properties
        sampleRate=100;   % how many points to record per second in time evolution
    end
    
    properties(Dependent)  % All link to coeffMatrix
        mesh;
        meshPars;
        simuPars;
        problemPars;
        M; S; MM; SS; MP; MG; vecQ; CB; vecR; MAB; vecF;  % Coeff matrices
        baseFunction;
    end
    
    methods
        function obj=FEMsolver3D(varargin)
            obj.coeffMatrix=PumpingDiffusionFEMSolver.CoeffMatrixGenerator(varargin{:});
        end
        
    end
    
    methods
        function plotMesh(obj,displayText)
            if nargin<2
                displayText=false;
            end
            obj.mesh.plotMesh(displayText);
        end
        function makeMesh(obj,forceRegenerate)
            if nargin<2
                forceRegenerate=false;
            end
            obj.mesh.makeMesh(forceRegenerate);
        end
        function genCoeffs(obj)
            obj.coeffMatrix.genCoeffs();
        end
        function cleanCoeffMatrices(obj)
            obj.coeffMatrix.cleanCoeffMatrices();
        end
        
        % Time evolution
        timeEvolution( obj,filename, u0 );
        % calc the expansion coeffs of initial state under the base function
        [ u0 ] = getInitialState(obj);
        

        
        % Save simulation result and parameters to file.
        [ filename ] = saveResultToFile( obj, filename );

    end
    
 
% ========================= Output Function for Time Evolution ===========================
    
    properties(Access=private)
        lastOutputTime=uint64(0);
        lastReportedT=0;
        lastSaveTime=uint64(0);
        lastSavedT=0;
    end
    properties(SetAccess=private)
        tyRecord=[];
    end
    
    methods(Access=private)

        % The Output Function for ode solvers
        % Output evolution progress and save result during evolution
        [ status ]=timeEvolutionOutputFunction(obj,t,y,flag);
        
    end
    
    
    
% ============================== Setters/Getters ========================================
    methods
% ------------------------- Setter Functions ----------------------------------------
        function set.meshPars(obj,meshPars)
            obj.coeffMatrix.mesh.meshPars=meshPars;
        end
        function set.simuPars(obj,simuPars)
            obj.coeffMatrix.simuPars=simuPars;
        end
        function set.problemPars(obj,problemPars)
            obj.coeffMatrix.problemPars=problemPars;
        end
% ------------------------- Getter Functions ----------------------------------------
        function mesh=get.mesh(obj)
            mesh=obj.coeffMatrix.mesh;
        end
        function meshPars=get.meshPars(obj)
            meshPars=obj.mesh.meshPars;
        end
        function simuPars=get.simuPars(obj)
            simuPars=obj.coeffMatrix.simuPars;
        end
        function problemPars=get.problemPars(obj)
            problemPars=obj.coeffMatrix.problemPars;
        end
        function baseFunction=get.baseFunction(obj)
            baseFunction=obj.coeffMatrix.baseFunHandle;
        end
        function val=get.M(obj)
            val=obj.coeffMatrix.M;
        end
        function val=get.MM(obj)
            val=obj.coeffMatrix.MM;
        end
        function val=get.S(obj)
            val=obj.coeffMatrix.S;
        end
        function val=get.SS(obj)
            val=obj.coeffMatrix.SS;
        end
        function val=get.MP(obj)
            val=obj.coeffMatrix.MP;
        end
        function val=get.MG(obj)
            val=obj.coeffMatrix.MG;
        end
        function val=get.MAB(obj)
            val=obj.coeffMatrix.MAB;
        end
        function val=get.CB(obj)
            val=obj.coeffMatrix.CB;
        end
        function val=get.vecQ(obj)
            val=obj.coeffMatrix.vecQ;
        end
        function val=get.vecR(obj)
            val=obj.coeffMatrix.vecR;
        end
        function val=get.vecF(obj)
            val=obj.coeffMatrix.vecF;
        end
    end
    
end

