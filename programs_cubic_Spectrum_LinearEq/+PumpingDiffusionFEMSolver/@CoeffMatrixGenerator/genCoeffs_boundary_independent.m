function [  ] = genCoeffs_boundary_independent( obj )
%genCoeffs_boundary_independent  Generate MM, SS, MG matrix as well as M, S, MP

    % get mesh info
    mesh=obj.mesh;
    Ndomains=mesh.Ndomains;
    K=obj.simuPars.maxOrder;
    baseFun=obj.baseFunHandle;
    Nbasis=baseFun.Nbasis;
    getNoByIxyz=baseFun.getNoByIxyz;  % getNoByIxyz should not contain zeros. i.e. all base should be used.
                                      % otherwise some reference below such as
                                      %    M(getNoByIxyz(1:K,1:K,1:K,Did),getNoByIxyz(1:K,1:K,1:K,Did))
                                      % will return errors.

    % ---------------------------- Generate M and S --------------------------------------------

    %Ixyz2No=reshape(1:K^3,K,K,K);     % basis mapping: (ix,iy,iz) -> rowNo/colNo in 
                                      %              mmm, mms, msm, smm and MP_kron

    nnzM=nnz(baseFun.phiphi(1:K,1:K,1));
    nnzS=nnz(abs(baseFun.dphidphi(1:K,1:K,1))+abs(baseFun.phiphi(1:K,1:K,1)));
                                      
    M=spalloc(Nbasis,Nbasis,Ndomains*(nnzM+2)^3+1);
    S=spalloc(Nbasis,Nbasis,Ndomains*(nnzS+2)^3+1);
    % enummerate domains
    for Did=1:Ndomains
        h=mesh.domains.h(1:3,Did);
        % get inner product between base
        phiphi_x  =baseFun.phiphi  (1:K,1:K,h(1));    % these are sparse matrices
        dphidphi_x=baseFun.dphidphi(1:K,1:K,h(1));
        phiphi_y  =baseFun.phiphi  (1:K,1:K,h(2));
        dphidphi_y=baseFun.dphidphi(1:K,1:K,h(2));
        phiphi_z  =baseFun.phiphi  (1:K,1:K,h(3));
        dphidphi_z=baseFun.dphidphi(1:K,1:K,h(3));
        % kron them
        mmm=kron(kron(phiphi_z  ,phiphi_y  ),phiphi_x  );
        mms=kron(kron(dphidphi_z,phiphi_y  ),phiphi_x  );
        msm=kron(kron(phiphi_z  ,dphidphi_y),phiphi_x  );
        smm=kron(kron(phiphi_z  ,phiphi_y  ),dphidphi_x);
        
        id=reshape(getNoByIxyz(1:K,1:K,1:K,Did),K^3,1);
        
        M(id,id)=M(id,id)+mmm; %#ok<SPRIX>
          
        S(id,id)=S(id,id)+mms+msm+smm; %#ok<SPRIX>
    end

    % -------------------------- Generate MP ---------------------------------------------------
    % assume that the weight function P(r) is separable as P(r)=A(x)B(y)C(z)
    % the MP_kron below is of the same size of MP when Ndomains=1, thus may require a lot more memory!
    % calc weighted inner product
    weightFun=obj.problemPars.funP;
    
    Nd_threshold=10;   % if Ndomains<Nd_threshold, use full matrix for MP
        
    if Ndomains<Nd_threshold
        MP=zeros(Nbasis,Nbasis);
    else
        MP=spalloc(Nbasis,Nbasis,round(0.02*Nbasis*Nbasis));  % The estimation of non-zero elements is not accurate
    end
    for Did=1:Ndomains
        xyz=mesh.domains.xyz(:,1,Did);
        h=mesh.domains.h(:,Did);
        % use kron product to generate MP
        phiphiF_x=h(1)/2*baseFun.innerProduct(@(x)weightFun{1}( (x+1)/2*h(1)+xyz(1) )...
                                        ,K,2*K+obj.simuPars.ngp);
        phiphiF_y=h(2)/2*baseFun.innerProduct(@(y)weightFun{2}( (y+1)/2*h(2)+xyz(2) )...
                                        ,K,2*K+obj.simuPars.ngp);
        phiphiF_z=h(3)/2*baseFun.innerProduct(@(z)weightFun{3}( (z+1)/2*h(3)+xyz(3) )...
                                        ,K,2*K+obj.simuPars.ngp);
        MP_kron=kron(kron(phiphiF_z,phiphiF_y),phiphiF_x);

        id=reshape(getNoByIxyz(1:K,1:K,1:K,Did),K^3,1);
        if issparse(MP)
            MP_kron(abs(MP_kron)<1e-23)=0;   % Assume that max(abs(MP)) is much smaller than 1e+6
            MP_kron=sparse(MP_kron);
        end
        MP(id,id)=MP(id,id)+MP_kron;
    end
    clear MP_kron
    
    if max(abs(MP(:)))>1e5  % check assumption
        error('The code here assume that max(abs(MP)) is much smaller than 1e+6.');
    end
    
    eps_MP=max(abs(MP(:)))*1e-17;
    disp(['All matrix elements of MP that smaller than ',num2str(eps_MP),' will be set to zero!']);
    if issparse(MP)
        [rowI,colI,val]=find(MP);
        val(abs(val)<eps_MP)=0;
        MP=sparse(rowI,colI,val);
    else
        MP(abs(MP)<eps_MP)=0;
    end
    
    
    % record results
    obj.MP=MP;
    obj.M=M;
    obj.S=S;
end

