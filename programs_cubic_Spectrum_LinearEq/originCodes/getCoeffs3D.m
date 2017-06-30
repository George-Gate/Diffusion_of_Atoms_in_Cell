function [ M, MM, SS, MG, Nbasis, fun2No, No2fun, getNoByIxyz, vecQ ] = getCoeffs3D( basis, mesh, K, dimRho, P, C0, D0, boundaryType )
% Generate Coefficient matrices for FEM
% [Inputs]
% basis: 'Bilinear+Lobatto'
% mesh: cuboid mesh generated by makeMesh3D()
% K   : cutoff of the order of Lobatto polynomials
% dimRho: the dimension of rho
% P   : 3x1 cell array, each element is a function handle for generate matrix MG
%       P{1}=@(x)...; P{2}=@(y)...; P{3}=@(z)...;   should have relation P(r)=P{1}*P{2}*P{3};
% C0, D0 : (dimRho x dimRho) matrix, see documentation for details
% boundaryType: can be 'first', 'firstUniform', 'second', 'secondZero'
%
% [Outputs]
% MM, SS, MG, vecQ: matrices for evolution. See documentation
% 
%
use_mex=0;
    switch basis
        case 'Bilinear+Lobatto'
            if use_mex
                [Mvec, Svec, fun2IntegralID_MP, MPid, Nbasis, fun2No, No2fun, getNoByIxyz]=getCoeffs3D_Lobatto_mex( boundaryType, mesh, K );
            else
                [Mvec, Svec, fun2IntegralID_MP, MPid, Nbasis, fun2No, No2fun, getNoByIxyz]=getCoeffs3D_Lobatto( boundaryType, mesh, K );
            end
            M=sparse(Mvec(:,1),Mvec(:,2),Mvec(:,3),Nbasis,Nbasis);
            S=sparse(Svec(:,1),Svec(:,2),Svec(:,3),Nbasis,Nbasis);
            MM=M; SS=S;
            for i=1:dimRho-1
                MM=blkdiag(MM,M);
                SS=blkdiag(SS,S);
            end
            % --------------- calc MP --------------------------
            ngp=2*K+100;
            MP=calcMP(mesh,P,fun2IntegralID_MP,MPid,K,Nbasis,ngp);
            
            matStat(MM);matStat(SS);matStat(MP);
            
            % --------------- calc MG --------------------------
            matDensity=nnz(MP)/numel(MP);
            if (matDensity*0.5938>0.2)
                % use full matrix to save memory
                MG=kron(full(C0),full(MP));
                MG_2=kron(sparse(D0),M);
                id=find(MG_2);
                MG(id)=MG(id)+MG_2(id);
                clear id MG_2
            else
                MG=kron(sparse(C0),MP)+kron(sparse(D0),M);
            end
            matStat(MG);
            
            % --------------- calc vecQ ------------------------
            vecQ=zeros(Nbasis*dimRho,1);
            if strcmp(boundaryType,'second')
                % calc vecQ
            end
            
        otherwise
            error(['Unknow basis: ',basis]);
    end


end