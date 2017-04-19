function [ Mvec, Svec, MPintegralList, QintegralList  ] = getCoeff_Lobatto( boundaryType, mesh )
% Generate the coeff matrices for a box area.
%    Do not support non constant first kind boundary condition so far. 
%    mesh.xList; mesh.yList; mesh.zList :   should be col vectors
%    mesh.type='cubic'
%    boundaryType: 'first', 'second', 'firstConstant', 'secondZero'
%
%
%
%   MPintegralList: a list of integral to be calculated. matrix MP depends on these integrals.
%     MPintegralList(i,k).ndim
%                        .ub/.lb
%                        .fun :     @(weightFun,x,y,z,...)
%      We suppose there is a function that can calc the integral: 
%             int(@(x,y,z...).fun(weightFun,x,y,z,...),lb,ub);
%      Where weightFun is a ndim-D function given by user
%   
%   QintegralList: similar to MPintegralList, its value give vecQ
%
%
%
%


% read mesh
if strcmp(boundaryType,'firstConstant')
    boundaryBasis=0;
else
    boundaryBasis=1;
end
switch mesh.type
    case 'cubic'
        Nx=length(mesh.xList);
        Ny=length(mesh.yList);
        Nz=length(mesh.zList);
        if (boundaryBasis)
            Nnodes=(Nx)*(Ny)*(Nz);
            Nedges=
            Nsurfaces=
        else
            Nnodes=(Nx-2)*(Ny-2)*(Nz-2);
            Nedges=
            Nsurfaces=
        end
        Ndomains=(Nx-1)*(Ny-1)*(Nz-1);
        
        
    otherwise
        error(['Unknow mesh type: ',mesh.type]);
end
numel(MPintegralList);
end

