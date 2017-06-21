classdef SimuPars
    %SimuPars   Provide simulation parameters for @FEMsolver3D, @@CoeffMatrixGenerator
    
    properties
        maxOrder=5;               % the maximun order of base on each direction
        baseFunHandle      % the base function object
        ngp=100;            % default number of gauss points when calculating numerical integration
    end
    
    properties(Dependent)
        basisName;       % the name of basis, read from baseFunHandle
    end
    
    methods
        function obj=SimuPars(baseFunHandle)
            if (nargin==0)
                obj.baseFunHandle=PumpingDiffusionFEMSolver.LobattoBase();
            elseif (nargin==1)
                obj.baseFunHandle=baseFunHandle;
            else
                error('Invalid number of input arguments.');
            end
        end
        
        function yes=eq(o1,o2)
            yes=true;
            if ~isa(o1,'PumpingDiffusionFEMSolver.SimuPars') || ...
                    ~isa(o2,'PumpingDiffusionFEMSolver.SimuPars')
                yes=false;
                return;
            end
            if ~(o1.maxOrder==o2.maxOrder && o1.ngp==o2.ngp && strcmp(o1.basisName,o2.basisName))
                yes=false;
                return;
            end
            if strcmp(o1.basisName,'')
                yes=false;
                return;
            end
        end
        
        function obj=set.baseFunHandle(obj,val)
            if isa(val,'PumpingDiffusionFEMSolver.BaseFunction')
                obj.baseFunHandle=val;
            else
                error('baseFunHandle should be a @PumpingDiffusionFEMSolver.BaseFunction');
            end
        end
        function basisName=get.basisName(obj)
            if isa(obj.baseFunHandle,'PumpingDiffusionFEMSolver.BaseFunction')
                basisName=obj.baseFunHandle.basisName;
            else
                basisName='';
            end
        end
    end
end

