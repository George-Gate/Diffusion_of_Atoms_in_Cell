function [  ] = genVecQ( obj )
%genVecQ  Generate vecQ for sceond type boundary condition
    
    obj.vecQ=obj.baseProjectionOnBoundary(obj.problemPars.drho_b,...
                         @(x,i)obj.baseFunHandle.funVal(x,i),   0);
    % multiply the diffusion constant
    obj.vecQ=obj.vecQ*obj.problemPars.D;
    
end

