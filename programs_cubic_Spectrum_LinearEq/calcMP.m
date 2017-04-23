function [MP, gp_x, gw, phi_num] = calcMP( mesh,funP,fun2IntegralID_MP,MPid,K,Nbasis,ngp )
% Calculate MP matrix for a given P function 
% funP should be separable: P(r)=A(x)*B(y)*C(z)
% Use Gaussian quadrature
% [Inputs]
% ngp: number of gauss points to be used for integration
% funP: cell array of function handle, P{1:3} for A(x), B(y) and C(z)

% --------------- Generate gauss points and weights --------------
% gp_x: gauss points; gw: weights
    [ gp_x, gw ] = getGaussPts( ngp );
% ----------------- Evaluate basis ----------------------------
    phi_num=zeros(ngp,K+2);
    phi_num(:,1)=(1-gp_x)/2;
    phi_num(:,2)=(1+gp_x)/2;
    phi_num(:,3:K+2)=lobattoP_N(1:K,gp_x);
% ---------------- Calc integrals -----------------------------
    IntegralList_MP=zeros(max(fun2IntegralID_MP(:)),1);
    for Did=1:mesh.Ndomains
        for dir=1:3
            intRange=mesh.domains.xyz(dir,:,Did)';
            P_num=gw.*funP{dir}( (gp_x+1)/2*(intRange(2)-intRange(1))+intRange(1) );
            P_num=repmat(P_num,1,K+2);
            tmpIntGrid=( intRange(2)-intRange(1) )/2*phi_num'*(P_num.*phi_num); % tmpIntGrid(l+1,m+1)=(phi_l,funP*phi_m);
            IntegralList_MP(fun2IntegralID_MP(:,:,dir,Did))=tmpIntGrid(:,:);
        end
    end
% -------------- Combine MP -----------------------------------
    MPvec=IntegralList_MP(MPid(:,3)).*IntegralList_MP(MPid(:,4)).*IntegralList_MP(MPid(:,5));
    if (mesh.Ndomains<63)
        MP=vec2full(MPid(:,1),MPid(:,2),MPvec,Nbasis,Nbasis);
    else
        MP=sparse(MPid(:,1),MPid(:,2),MPvec,Nbasis,Nbasis);
    end
    if (mesh.Ndomains>5)
        MP=sparse(MP);
    end
%     toc;tic;
%     MP2=sparse(MPid(:,1),MPid(:,2),MPvec,Nbasis,Nbasis);
%     toc;
%     max(max(abs(MP2-MP)))
end

