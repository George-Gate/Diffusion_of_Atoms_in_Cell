classdef Mesh3D < handle
    %Mesh3D  Generate a square mesh for a cubic domain
    %   Detailed explanation goes here
    
    properties
        use_mex=1;   % indicate whether to use mex function for makeMesh() method or not.
        meshPars;    % parameters for mesh, @MeshPars
    end
    
    properties(SetAccess=private)
        meshPars_current=0;  % the meshPars for the currently generated mesh
        nodes;
        edges;
        surfaces;
        domains;
        Nnodes;
        Nedges;
        Nsurfaces;
        Ndomains;
    end
    
    properties(Dependent, SetAccess=private)
        outDated;     % indicate whether the current mesh is generated from current meshPars.
    end

    
    methods
        % constructor
        function obj=Mesh3D(meshPars)
            import PumpingDiffusionFEMSolver.MeshPars
            if nargin==0
                obj.meshPars=MeshPars();
            else
                if isa(meshPars,'MeshPars')
                    obj.meshPars=meshPars;
                else
                    error('meshPars should be an instance of class PumpingDiffusionFEMSolver.MeshPars');
                end
            end
        end
        
        % setter & getter functions
        function set.meshPars(obj,val)
            if isa(val,'PumpingDiffusionFEMSolver.MeshPars')
                obj.meshPars=val;
            else
                error('Input should be an object of class @PumpingDiffusionFEMSolver.MeshPars');
            end
        end
        function val=get.outDated(obj)
            val= ~(obj.meshPars_current == obj.meshPars );
        end
        
        
        % generate mesh
        function makeMesh(obj,forceRegenerate)
            if nargin<2
                forceRegenerate=false;
            end
            if forceRegenerate || obj.outDated
                startT=tic;
                import PumpingDiffusionFEMSolver.Mesh3D
                % call function makeMesh3D_cubic() to get mesh structure
                if obj.use_mex
                    mesh0=Mesh3D.makeMesh3D_cubic_mex(obj.meshPars.xList,  obj.meshPars.yList,  obj.meshPars.zList);
                else
                    mesh0=Mesh3D.makeMesh3D_cubic(obj.meshPars.xList,  obj.meshPars.yList,  obj.meshPars.zList);
                end
                % save mesh structure
                obj.nodes=mesh0.nodes;
                obj.edges=mesh0.edges;
                obj.surfaces=mesh0.surfaces;
                obj.domains=mesh0.domains;
                obj.Nnodes=mesh0.Nnodes;
                obj.Nedges=mesh0.Nedges;
                obj.Nsurfaces=mesh0.Nsurfaces;
                obj.Ndomains=mesh0.Ndomains;
                % set meshPars_current
                obj.meshPars_current=obj.meshPars;
                disp(['Time used to generate mesh: ',sec2hms(toc(startT))]);
            else
                % meshPars not changed, skip
                disp('meshPars wasn''t changed since last makeMesh. Use .makeMesh(1) to force regenerate.');
            end
        end
        
        % Visualize mesh
        % displayText and displayRelation is two logical flags that are optional
        plotMesh(obj,displayText);

    end
    
    methods(Hidden, Access=private, Static)
        % the actual function that generate a mesh for cubic domain. Can be complied to an mex-function. 
        mesh = makeMesh3D_cubic(xList,yList,zList);
        mesh = makeMesh3D_cubic_mex(xList,yList,zList);
    end
    
end

