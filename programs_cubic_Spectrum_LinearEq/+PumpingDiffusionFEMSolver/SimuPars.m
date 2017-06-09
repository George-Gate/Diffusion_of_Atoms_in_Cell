classdef SimuPars
    %SimuPars   Provide simulation parameters for @FEMsolver3D
    %   Detailed explanation goes here
    
    properties
        K=5;               % the maximun order of base on each direction
        basisName
        baseFunHandle
        ngp=10;            % default number of gauss points when calculating numerical integration
    end
    
    methods
        function obj=SimuPars(basisName,baseFunHandle)
            if (nargin==0)
                obj.basisName='Lobatto';
                obj.baseFunHandle=PumpingDiffusionFEMSolver.Lobatto();
            elseif isa(baseFunHandle,PumpingDiffusionFEMSolver.BaseFunction)
                obj.basisName=basisName;
                obj.baseFunHandle=baseFunHandle;
            else
                error('Invalid number of input arguments or the type of baseFunHandle is invalid.');
            end
        end
        function obj=set.baseFunHandle(obj,val)
            if isa(val,PumpingDiffusionFEMSolver.BaseFunction)
                obj.baseFunHandle=val;
            else
                error('baseFunHandle should be a @PumpingDiffusionFEMSolver.BaseFunction.');
            end
        end
    end
    
end

