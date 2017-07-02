function [ filename ] = saveResultToFile( obj, filename )
%saveResultToFile  Save simulation result and parameters to file
% Also set obj.FEMResult
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
end

