function Gf = calcGf( mesh,funP,constf,matC,D0,getNoByIxyz,dimRho,K,Nbasis,ngp )
% Calculate vector Gf for a given P function and matrix C0, D0
% funP should be separable: P(r)=A(x)*B(y)*C(z)
% Use Gaussian quadrature
% See the 1.1 section (constant first kind boundary condition) of algorithm description.
%
%      Gf -> (v_i, funP)*(matC*constf) + (v_i, 1)*(matD*constf)
%
% [Inputs]
% ngp: number of gauss points to be used for integration
% funP: cell array of function handle, P{1:3} for A(x), B(y) and C(z)
% constf: a column vector, the boundary condition
% D0: alpha*T0*matD

% ---------------- Calc matD part integrals -----------------------------
    
    Gf=zeros(dimRho*Nbasis,1);
    D0_f=D0*constf;
    vec1=basisInnerProduct( 1, getNoByIxyz, K, Nbasis, mesh, ngp);
    for j=1:dimRho
        Gf((j-1)*Nbasis+1:j*Nbasis)=D0_f(j)*vec1;
    end

% ---------------- Calc matC part integrals -----------------------------
    % vecP(No)=(v_No,P(r))
    vecP=basisInnerProduct( funP, getNoByIxyz, K, Nbasis, mesh, ngp);
    C0_f=matC*constf;
    Gf=Gf+kron(C0_f,vecP);
end

