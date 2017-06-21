function [  ] = genMABvecF( obj )
%genMABvecF   Generate MAB and vecF for robin type boundary condition
    
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
    ngp=2*maxOrder+obj.simuPars.ngp;
% -------------------- Calc MAB ---------------------------------------------
    A=obj.problemPars.robinA;
    obj.MAB=spalloc(Nbasis*obj.dimRho,Nbasis*obj.dimRho,Nbasis*maxOrder*100);
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
                
                % get matrix A for the current surface
                if isnumeric(A)
                    A0=A;
                else
                    A0=A{dkDir2lnID(dk,dir)};
                end
                
% =============== if A0 is constant or has a matrix-function separable form =================
                % then use direct product to speed up
                if isnumeric(A0) || numel(A0)<=2
               % ---------- Set A1: MAB=kron(A1,MA0) ------------
                    if isscalar(A0)  % A0 is a 2D function or a constant scalar
                        A1=eye(obj.dimRho);
                    elseif isnumeric(A0)  % A0 is a constant matrix
                        A1=A0;
                    elseif isnumeric(A0{1})  % A0 has the matrix-function separable form
                        A1=A0{1};        % then A0{1} should be a constant matrix
                    else
                        error('The format of boundary condition is invalid.');
                    end

               % ---------- calc base value on boundary ---------------
                    int3=base.funVal(planePos,1:maxOrder);
                    int3=sparse(int3.*int3.');
               % ---------- calc weighted inner product int12 ---------
                    % if A0 contains a 2D function
                    if iscell(A0) && (isa(A0{1},'function_handle') || (length(A0)==2 && isa(A0{2},'function_handle')))
                        if isa(A0{1},'function_handle')  % A0 is a scalar separable 2D function
                            funA=A0{1};
                        else        % A0 has the matrix-function separable form
                            funA=A0{2};
                        end
                        if ~exist('gp_y','var') || length(gp_y)~=ngp   % reuse gp_y and gw for speed up
                            import PumpingDiffusionFEMSolver.library.getGaussPts
                            [ gp_y, gw ] = getGaussPts( ngp );
                        end
                        int1=zeros(maxOrder,maxOrder,ngp);
                        % calc integral along axis 1
                        for iP=1:ngp
                            y=(gp_y(iP)+1)*h2/2+range2(1);
                            int1(:,:,iP)=gw(iP)*h1/2*base.innerProduct( @(x)funA( (x+1)*h1/2+range1(1) , y ) ,...
                                                  maxOrder,ngp);
                        end
                        % calc base fun value at gauss points
                        funVal2=base.funVal(gp_y,1:maxOrder);
                        funVal2=permute(funVal2,[2,3,1]).*permute(funVal2,[3,2,1]);
                        % calc 2-D integral
                        int12=zeros(maxOrder^2,maxOrder^2,ngp);
                        for iP=1:ngp
                            int12(:,:,iP)=kron(funVal2(:,:,iP),int1(:,:,iP));
                        end
                        int12=h2/2*sum(int12,3);    % int12((i1,i2),(k1,k2)) -> int( v_i1(x)*v_i2(y)*v_k1(x)*v_k2(y)*f(x,y) )

                    % if A0 is constant or contains separable 2D function    
                    else
                        % set funA
                        funA=cell(2,1);
                        if isnumeric(A0) && isscalar(A0)  % A0 constant
                            funA{1}=A0;funA{2}=1;
                        elseif isnumeric(A0)            % A0 constant matrix
                            funA{1}=1;funA{2}=1;
                        else  % A0 is cell
                            if isnumeric(A0{1})   % A0 matrix-function separable form
                                funA{1}=A0{2}{1};funA{2}=A0{2}{2};
                            else                  % scalar separable 2D function
                                funA{1}=A0{1}{1};funA{2}=A0{1}{2};
                            end
                        end
                        % calc int12
                        K=maxOrder;
                        if isnumeric(funA{1})
                            int1=funA{1}*base.phiphi(1:K,1:K,h1);
                            int2=funA{2}*base.phiphi(1:K,1:K,h2);
                        else
                            int1=h1/2*base.innerProduct(@(x)funA{1}((x+1)*h1/2+range1(1)),K,ngp);
                            int2=h2/2*base.innerProduct(@(y)funA{2}((y+1)*h2/2+range2(1)),K,ngp);
                        end
                        int12=kron(int2,int1);
                    end
               % ---------- calc MA0 ----------------
                    int0=kron(int3,int12);
                    MA0=spalloc(Nbasis,Nbasis,nnz(int0));
                    K=maxOrder;
                    getNoByIxyz_new=permute(base.getNoByIxyz(1:K,1:K,1:K,Did),[axis1(dir),axis2(dir),axis3(dir)]);
                    id=reshape(getNoByIxyz_new,K^3,1);
                    MA0(id,id)=int0; %#ok<SPRIX>
               % ---------- add to obj.MAB ----------------
                    obj.MAB=obj.MAB+kron(A1,MA0);
                    
                    clear getNoByIxyz_new MA0 int0 int12 int3 funVal2
                    
% =============== if A0 does not have the form f(x,y,z)*A1 ============================
                % then calc MAB element-wise
                elseif iscell(A0) && size(A0,1)==obj.dimRho && size(A0,2)==obj.dimRho
               % ---------- calc base value on boundary ---------------
                    int3=base.funVal(planePos,1:maxOrder);
                    int3=sparse(int3.*int3.');
                    K=maxOrder;
                    getNoByIxyz_new=permute(base.getNoByIxyz(1:K,1:K,1:K,Did),[axis1(dir),axis2(dir),axis3(dir)]);
                    id=reshape(getNoByIxyz_new,K^3,1);
               % ---------- loop over all elements of cell array A0 -------------
                    % A0{i,j} has three possibility: number, 2D function, cell array containing two 1D functions
                    for rID=1:obj.dimRho
                        for cID=1:obj.dimRho
                        % ------------ Calc int12 -------------------------
                            % A0{rID,cID} is number or separable 2D function
                            if isnumeric(A0{rID,cID}) || iscell(A0{rID,cID})
                                K=maxOrder;
                                if isnumeric(A0{rID,cID})
                                    int1=A0{rID,cID}*base.phiphi(1:K,1:K,h1);
                                    int2=base.phiphi(1:K,1:K,h2);
                                else
                                    int1=h1/2*base.innerProduct(@(x)A0{rID,cID}{1}((x+1)*h1/2+range1(1)),K,ngp);
                                    int2=h2/2*base.innerProduct(@(y)A0{rID,cID}{2}((y+1)*h2/2+range2(1)),K,ngp);
                                end
                                int12=kron(int2,int1);
                                
                            % A0{rID,cID} is 2D function
                            elseif isa(A0{rID,cID},'function_handle')
                                if ~exist('gp_y','var') || length(gp_y)~=ngp   % reuse gp_y and gw for speed up
                                    import PumpingDiffusionFEMSolver.library.getGaussPts
                                    [ gp_y, gw ] = getGaussPts( ngp );
                                end
                                int1=zeros(maxOrder,maxOrder,ngp);
                                % calc integral along axis 1
                                for iP=1:ngp
                                    y=(gp_y(iP)+1)*h2/2+range2(1);
                                    int1(:,:,iP)=gw(iP)*h1/2*base.innerProduct( @(x)A0{rID,cID}( (x+1)*h1/2+range1(1) , y ) ,...
                                                          maxOrder,ngp);
                                end
                                % calc base fun value at gauss points
                                funVal2=base.funVal(gp_y,1:maxOrder);
                                funVal2=permute(funVal2,[2,3,1]).*permute(funVal2,[3,2,1]);
                                % calc 2-D integral
                                int12=zeros(maxOrder^2,maxOrder^2,ngp);
                                for iP=1:ngp
                                    int12(:,:,iP)=kron(funVal2(:,:,iP),int1(:,:,iP));
                                end
                                int12=h2/2*sum(int12,3);    % int12((i1,i2),(k1,k2)) -> int( v_i1(x)*v_i2(y)*v_k1(x)*v_k2(y)*f(x,y) )

                            else
                                error('The format of boundary condition is invalid.');
                            end
                        % ------------ Calc MAij ----------------------------
                            int0=kron(int3,int12);
                            MAij=spalloc(Nbasis,Nbasis,nnz(int0));
                            MAij(id,id)=int0; %#ok<SPRIX>
                        % ------------ add to obj.MAB ------------------------
                            rL=(rID-1)*Nbasis+(1:Nbasis);
                            cL=(cID-1)*Nbasis+(1:Nbasis);
                            obj.MAB(rL,cL)=obj.MAB(rL,cL)+MAij;
                        end
                    end
                    clear rL cL id getNoByIxyz_new MAij int0 int12 int3 funVal2
                else
                    error('The format of boundary condition is invalid.');
                end
                
            end
        end
    end
    
    % delete the pre-allocated space for obj.MAB
    MAB_old=obj.MAB;
    obj.MAB=spalloc(size(MAB_old,1),size(MAB_old,2),nnz(MAB_old));
    obj.MAB(:)=MAB_old(:);

% -------------------- Calc vecF ---------------------------------------------
    obj.vecF=obj.baseProjectionOnBoundary(obj.problemPars.robinF,...
                      @(x,i)obj.baseFunHandle.funVal(x,i),   0);
    % multiply the diffusion constant
    obj.vecF=obj.vecF*obj.problemPars.D;


end

