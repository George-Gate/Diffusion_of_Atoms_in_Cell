function [  ] = genCoeffs_boundary_independent_slow( obj )
%genCoeffs_boundary_independent  Generate MM, SS, MG matrix as well as M, S, MP

    % get mesh info
    mesh=obj.mesh;
    Ndomains=mesh.Ndomains;
    K=obj.simuPars.maxOrder;
    baseFun=obj.baseFunHandle;
    getNoByIxyz=baseFun.getNoByIxyz;

    % ---------------------------- Generate Mvec and Svec --------------------------------------------
    % get inner product between base
    phiphi_h=full(baseFun.phiphi(1:K,1:K,1));    % these two are sparse matrix
    dphidphiH=full(baseFun.dphidphi(1:K,1:K,1)); 

    Mvec=zeros(Ndomains*(nnz(phiphi_h )+2)^3+1,3);
    Svec=zeros(Ndomains*(nnz(dphidphiH)+2)^3+1,3);
    topM=1; topS=1;

    % enummerate basis to get Mvec and Svec
    for Did=1:Ndomains
        h=mesh.domains.h(1:3,Did);
        [iList,kList]=find(phiphi_h~=0 | dphidphiH~=0);
        for xID=1:length(iList)   % [row, column] for x basis
            ix=iList(xID);
            kx=kList(xID);
    %       disp([num2str([ix,kx])]);
            phiphi_x=phiphi_h(ix,kx)*h(1);
            dphidphi_x=dphidphiH(ix,kx)/h(1);
            for yID=1:length(iList)  % [row, column] for y basis
                iy=iList(yID);
                ky=kList(yID);

                phiphi_y=phiphi_h(iy,ky)*h(2);
                dphidphi_y=dphidphiH(iy,ky)/h(2);
                for zID=1:length(iList)    % [row, column] for z basis
                    iz=iList(zID);
                    kz=kList(zID);
                    phiphi_z=phiphi_h(iz,kz)*h(3);
                    dphidphi_z=dphidphiH(iz,kz)/h(3);

                    rowNo=getNoByIxyz(ix,iy,iz,Did);
                    colNo=getNoByIxyz(kx,ky,kz,Did);

                    % set M matrix
                    Mvec(topM,3)=phiphi_x*phiphi_y*phiphi_z;
                    if abs(Mvec(topM,3))>0
                        Mvec(topM,1)=rowNo;
                        Mvec(topM,2)=colNo;
                        topM=topM+1;
                    end
                    % set S matrix
                    Svec(topS,3)=dphidphi_x*phiphi_y*phiphi_z ...
                                +phiphi_x*dphidphi_y*phiphi_z ...
                                +phiphi_x*phiphi_y*dphidphi_z;
                    if abs(Svec(topS,3))>0
                        Svec(topS,1)=rowNo;
                        Svec(topS,2)=colNo;
                        topS=topS+1;
                    end
                end
            end
        end
    end

    % trim Mvec & Svec
    Mvec=Mvec(1:topM-1,1:3);
    Svec=Svec(1:topS-1,1:3);
    Nbasis=baseFun.Nbasis;
    M=sparse(Mvec(:,1),Mvec(:,2),Mvec(:,3),Nbasis,Nbasis);
    S=sparse(Svec(:,1),Svec(:,2),Svec(:,3),Nbasis,Nbasis);

    % -------------------------- Generate MP ---------------------------------------------------
    % assume that the weight function P(r) is separable as P(r)=A(x)B(y)C(z)
    % calc weighted inner product
    weightFun=obj.problemPars.funP;

    Ixyz2No=reshape(1:K^3,K,K,K);     % basis mapping: (ix,iy,iz) -> rowNo in MP_kron

    MP=zeros(K^3,K^3);
    for Did=1:Ndomains
        xyz=mesh.domains.xyz(:,1,Did);
        h=mesh.domains.h(:,Did);
        % use kron product to generate MP
        phiphiF_x=baseFun.innerProduct(@(x)weightFun{1}( (x+1)/2*h(1)+xyz(1) )...
                                        ,K,obj.simuPars.ngp);
        phiphiF_y=baseFun.innerProduct(@(y)weightFun{2}( (y+1)/2*h(2)+xyz(2) )...
                                        ,K,obj.simuPars.ngp);
        phiphiF_z=baseFun.innerProduct(@(z)weightFun{3}( (z+1)/2*h(3)+xyz(3) )...
                                        ,K,obj.simuPars.ngp);
        MP_kron=kron(kron(phiphiF_z,phiphiF_y),phiphiF_x);

        MP(getNoByIxyz(1:K,1:K,1:K,Did),getNoByIxyz(1:K,1:K,1:K,Did))...
             =MP(getNoByIxyz(1:K,1:K,1:K,Did),getNoByIxyz(1:K,1:K,1:K,Did))...
             +MP_kron(Ixyz2No(1:K,1:K,1:K),Ixyz2No(1:K,1:K,1:K));
    end

    eps_MP=max(abs(phiphiF_x(:)))*max(abs(phiphiF_y(:)))*max(abs(phiphiF_z(:)))*1e-17;
    disp(['All matrix elements of MP that smaller than ',num2str(eps_MP),' will be set to zero!']);
    MP(abs(MP)<eps_MP)=0;
    
    % record results
    obj.MP=MP;
    obj.M=M;
    obj.S=S;
end

