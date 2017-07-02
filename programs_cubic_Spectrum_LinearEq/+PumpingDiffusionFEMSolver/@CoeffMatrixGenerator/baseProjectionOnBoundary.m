function [ vecX ] = baseProjectionOnBoundary( obj, funX, baseFunEvaluator, baseDerivative  )
%baseProjectionOnBoundary  To calc the vecX's in boundary terms
%  For the inner product (v_i,funX) on boundary, set baseFunEvaluator=obj.baseFunHandle.funVal and baseDerivative=0
%  For the inner product (dv_i/dn,funX) on boundary, set baseFunEvaluator=obj.baseFunHandle.funFirstDerivative and baseDerivative=1
%  For the inner product (dv_i/dz,funX(x,y)) on boundary, set baseFunEvaluator=obj.baseFunHandle.funFirstDerivative and baseDerivative=0

% ------------ Constant Definitions -----------------------------------------
    dkDir2lnID=[1,3,5;     % dkDir2lnID(dk,dir)=index in funX.
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
    
    vecX=zeros(Nbasis*obj.dimRho,1);
            
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
                
                % get funX for the current surface
                if isnumeric(funX)
                    funX0=funX;
                else
                    funX0=funX{dkDir2lnID(dk,dir)};
                end
                
          % ==================== If funX0 is a constant vector ======================================
                % We only need to calc the integral once
                if isnumeric(funX0) && length(funX0)==obj.dimRho && iscolumn(funX0)
                    % calc the inner product of base and constant
                    % assume that BaseFunction is polynomials whose order are up to maxOrder.
                    prj{1}=h1/2*base.projection(@(x)1,maxOrder,maxOrder);
                    prj{2}=prj{1}/h1*h2;
                    % calc the value of base on boundary
                    prj{3}=reshape(baseFunEvaluator(planePos,1:maxOrder),maxOrder,1);
                    if baseDerivative==1
                        prj{3}=sign(planePos)*prj{3};
                    end
                    % contruct three axis base and add them to vecX
                    prj0=kron(kron(prj{zaxis(dir)},prj{yaxis(dir)}),prj{xaxis(dir)});
                    vecX_tmp=zeros(Nbasis,1);  K=maxOrder;
                    vecX_tmp(reshape(base.getNoByIxyz(1:K,1:K,1:K,Did),K^3,1))=prj0;
                    vecX=vecX+kron(funX0,vecX_tmp);
                    clear prj prj0 vecX_tmp
          % =================== If funX0 is cell array with function handles =======================
                elseif iscell(funX0) && length(funX0)==obj.dimRho  && iscolumn(funX0)
                    % calc vecX for each dimension of rho
                    for iDim=1:obj.dimRho
                        funX0_0=funX0{iDim};
                 % --------------- funX0{iDim} is a constant or separable 2-D function ------------------
                        if isnumeric(funX0_0) || iscell(funX0_0)
                            % calc the inner product of base and funX0_0
                            if isnumeric(funX0_0)
                                % assume that BaseFunction is polynomials whose order are up to maxOrder.
                                prj{1}=h1/2*base.projection(@(x)1,maxOrder,maxOrder);
                                prj{2}=prj{1}/h1*h2*funX0_0;
                            else
                                prj{1}=h1/2*base.projection(@(x)funX0_0{1}((x+1)*h1/2+range1(1)),maxOrder,ngp);
                                prj{2}=h2/2*base.projection(@(x)funX0_0{2}((x+1)*h2/2+range2(1)),maxOrder,ngp);
                            end
                            % calc the value of base on boundary
                            prj{3}=reshape(baseFunEvaluator(planePos,1:maxOrder),maxOrder,1);
                            if baseDerivative==1
                                prj{3}=sign(planePos)*prj{3};
                            end
                            % contruct three axis base and add them to vecX
                            prj0=kron(kron(prj{zaxis(dir)},prj{yaxis(dir)}),prj{xaxis(dir)});
                            K=maxOrder;
                            vecX((iDim-1)*Nbasis+reshape(base.getNoByIxyz(1:K,1:K,1:K,Did),K^3,1))...
                               =vecX((iDim-1)*Nbasis+reshape(base.getNoByIxyz(1:K,1:K,1:K,Did),K^3,1))...
                               +prj0;
                            clear prj prj0
                 % --------------- funX0{iDim} is a non-separable 2-D function ---------------------------
                        else
                            if ~exist('gp_y','var') || length(gp_y)~=ngp   % reuse gp_y and gw for speed up
                                import PumpingDiffusionFEMSolver.library.getGaussPts
                                [ gp_y, gw ] = getGaussPts( ngp );
                            end
                            prj1=zeros(maxOrder,ngp);
                            % calc integral along axis 1
                            for iP=1:ngp
                                y=(gp_y(iP)+1)*h2/2+range2(1);
                                prj1(:,iP)=gw(iP)*h1/2*base.projection( @(x)funX0_0( (x+1)*h1/2+range1(1) , y ) ,...
                                                     maxOrder,ngp);
                            end
                            % calc 2-D integral
                            funVal2=base.funVal(gp_y,1:maxOrder);
                            prj12=h2/2*prj1*funVal2;    % prj12(i1,i2) -> int( v_i1(x)*v_i2(y)*f(x,y) )
                            % contruct three axis base and add them to vecX
                            K=maxOrder;
                            % calc the base value on boundary
                            prj3=reshape(baseFunEvaluator(planePos,1:K),K,1);
                            if baseDerivative==1
                                prj3=sign(planePos)*prj3;
                            end
                            prj0=kron(prj3,reshape(prj12,K^2,1));   % this is something like prj0=kron(kron(prj3,prj2),prj1);
                            I2No=permute(base.getNoByIxyz(1:K,1:K,1:K,Did),[axis1(dir),axis2(dir),axis3(dir)]);
                            vecX( (iDim-1)*Nbasis+reshape(I2No(1:K,1:K,1:K),K^3,1) )...
                               =vecX( (iDim-1)*Nbasis+reshape(I2No(1:K,1:K,1:K),K^3,1) )...
                               +prj0;
                            clear prj1 funVal2 prj12 prj3 prj0 I2No
                        end
                    end
                  
          % =================== If funX0 is either a column vector or a column cell array =======================
                else
                    error('The format of boundary condition is invalid.');
                end
                
            end
        end
    end

end

