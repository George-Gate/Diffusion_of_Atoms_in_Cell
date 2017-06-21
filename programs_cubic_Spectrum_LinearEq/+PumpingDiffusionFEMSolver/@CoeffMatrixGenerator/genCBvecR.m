function [  ] = genCBvecR( obj )
%genCBvecR  Generate CB, vecR for first type boundary condition

    
% ------------ Constant Definitions -----------------------------------------
    dkDir2lnID=[1,3,5;     % dkDir2lnID(dk,dir)=index in drho_b.
                2,4,6];
    axis1=[2;3;1];      % axis1(dir)=the id for axis 1 if the plane direction is dir
    axis2=[3;1;2];      % axis2(dir)=the id for axis 2 if the plane direction is dir
    axis3=[1;2;3];      % axis3=dir, the plane is orthonal to axis 3
    xaxis=[3;2;1];      % xaxis(dir)=the axis number of x axis when plane direction is dir
    yaxis=[1;3;2];      
    zaxis=[2;1;3];
% ----------------------------------------------------------------------------
    mesh=obj.mesh;
    base=obj.baseFunHandle;
    Nbasis=base.Nbasis;
    maxOrder=obj.simuPars.maxOrder;
    ngp=maxOrder+obj.simuPars.ngp;
    
% -------------------- Calc matrix CB ----------------------------------------
    K=maxOrder;
    % optimized for lobatto basis, may need to change nzmax for other basis
    CB0=spalloc(Nbasis,Nbasis,mesh.Ndomains*K^3*K*16);   % alloc a large enough sparse matrix
    % enumerate all surface and find boundary surface.
    for Did=1:mesh.Ndomains
        for dir=1:3
            for dk=1:2
                % skip non-boundary syrfaces
                Sid=mesh.domains.s(dk,dir,Did);
                if ~mesh.surfaces.onB(Sid)
                    continue;
                end
                % find integral area by surfaces.n[1] and surfaces.n[3]
                n1=mesh.surfaces.n(1,Sid);
                n3=mesh.surfaces.n(3,Sid);
                xyz=[mesh.nodes.x([n1;n3]),mesh.nodes.y([n1;n3]),mesh.nodes.z([n1;n3])];
                range1=[xyz(2,axis1(dir));xyz(1,axis1(dir))];
                range2=[xyz(2,axis2(dir));xyz(1,axis2(dir))];
                planePos=xyz(1,dir);
                h1=range1(2)-range1(1);h2=range2(2)-range2(1);
                % calc phiphi matrix
                phiphi{1}=base.phiphi(1:K,1:K,h1);   % this is sparse matrix
                phiphi{2}=base.phiphi(1:K,1:K,h2);
                phiphi{3}=base.funVal(planePos,1:K).'  .*  base.funFirstDerivative(planePos,1:K) * sign(planePos);
                phiphi{3}=sparse(phiphi{3}+phiphi{3}.');
                phiphi0=kron(kron(phiphi{zaxis(dir)},phiphi{yaxis(dir)}),phiphi{xaxis(dir)});
                % add to CB_tmp
                id=reshape(base.getNoByIxyz(1:K,1:K,1:K,Did),K^3,1);
                CB0(id,id)=CB0(id,id)+phiphi0; %#ok<SPRIX>
            end
        end
    end
    % get obj.CB
    obj.CB=CB0;
    for i=1:obj.dimRho-1
        obj.CB=blkdiag(obj.CB,CB0);
    end
    
    
% -------------------- Calc vecR ---------------------------------------------
    obj.vecR=obj.baseProjectionOnBoundary(obj.problemPars.rho_b,...
                      @(x,i)obj.baseFunHandle.funFirstDerivative(x,i),  1);
   
    % multiply the diffusion constant
    obj.vecR=obj.vecR*obj.problemPars.D;


end

