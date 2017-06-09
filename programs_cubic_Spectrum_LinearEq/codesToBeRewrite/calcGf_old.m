function [Gf, gp_x, gw, phi_num] = calcGf_old( mesh,funP,constf,matC,D0,getNoByIxyz,dimRho,K,Nbasis,ngp )
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

% --------------- Generate gauss points and weights --------------
% gp_x: gauss points; gw: weights
    [ gp_x, gw ] = getGaussPts( ngp );
% ----------------- Evaluate basis ----------------------------
    phi_num=zeros(ngp,K+2);
    phi_num(:,1)=(1-gp_x)/2;
    phi_num(:,2)=(1+gp_x)/2;
    phi_num(:,3:K+2)=lobattoP_N(1:K,gp_x);
% ---------------- Calc matD part integrals -----------------------------
    Gf=zeros(dimRho*Nbasis,1);
    D0_f=D0*constf;
    for Did=1:mesh.Ndomains
        for ix=0:2
            for iy=0:2
                for iz=0:2
                    % (v_i, 1)=-5/6*(i-0.5).*(i-0.5)+29/24;  i=0/1: 1.0,  i=2: -2/3
                    No=getNoByIxyz(ix+1,iy+1,iz+1,Did);
                    if (No==0); continue; end
                    Iv=(-5/6*(ix-0.5)*(ix-0.5)+29/24)...
                      *(-5/6*(iy-0.5)*(iy-0.5)+29/24)...
                      *(-5/6*(iz-0.5)*(iz-0.5)+29/24);
                    for j=1:dimRho
                        Gf(No+(j-1)*Nbasis)=Gf(No+(j-1)*Nbasis)+Iv*D0_f(j);
                    end
                end
            end
        end
    end
% ---------------- Calc matC part integrals -----------------------------
    vPList=ones(Nbasis,1);   % vPList(No)=(v_No,P(r))
    for Did=1:mesh.Ndomains
        for dir=1:3
            intRange=mesh.domains.xyz(dir,:,Did)';
            P_num=gw.*funP{dir}( (gp_x+1)/2*(intRange(2)-intRange(1))+intRange(1) );
            tmpIntList=( intRange(2)-intRange(1) )/2*phi_num'*P_num; % tmpIntList(l+1)=(phi_l,funP);
            for l=0:K+1   % times tmpIntList into vPList
                if dir==1
                    NoList=getNoByIxyz(l+1,:,:,Did);
                elseif dir==2
                    NoList=getNoByIxyz(:,l+1,:,Did);
                else
                    NoList=getNoByIxyz(:,:,l+1,Did);
                end
                NoList=NoList(NoList>0);
                vPList(NoList)=vPList(NoList)*tmpIntList(l+1);
            end
        end
    end
    C0_f=matC*constf;
    Gf=Gf+kron(C0_f,vPList);
end

