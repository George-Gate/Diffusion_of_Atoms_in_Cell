function [ u0 ] = getInitialState(obj,rho_0)
%getInitialState  Calc the expansion coeffs of initial state under the base function
%  Should be called after calling genCoeffs() method. Otherwise the status of baseFunction
%   may be not up to date.
    Nbasis=obj.baseFunction.Nbasis;
    dimRho=obj.coeffMatrix.dimRho;
    if nargin<2
        rho_0=obj.problemPars.rho_0;  % get dependent property rho_0 may be time expensive, store it.
    end
    base=obj.baseFunction;
    maxOrder=obj.simuPars.maxOrder;
    ngp=maxOrder+obj.simuPars.ngp;
    
    % check the form of rho_0 (by using the checkInitValFormat() method of class @ProblemPars)
    rhoType=checkInitFormat(obj.problemPars,rho_0);
    if rhoType==0
        error('The format of rho_0 is invalid.');
    elseif rhoType==1
        rho_0={rho_0};
    end
    v0=zeros(Nbasis*dimRho,1);
    for Did=1:obj.mesh.Ndomains   % Did for Domain ID
        % get integrate area
        h=obj.mesh.domains.h(:,Did);
        xyz0=obj.mesh.domains.xyz(:,1,Did);
        
        K=maxOrder;
        id=reshape(base.getNoByIxyz(1:K,1:K,1:K,Did),K^3,1);
        
        for iPart=1:length(rho_0)   % rho_0 is the summation of several parts
            rho_0_i=rho_0{iPart};
            for iDim=1:dimRho
                if isnumeric(rho_0_i)
                    fun={@(x)rho_0_i(iDim),@(y)1, @(z)1};
                elseif isnumeric(rho_0_i{iDim})
                    fun={@(x)rho_0_i{iDim},@(y)1, @(z)1};
                else   % rho_0_i{iDim} is 3D function
                    fun=rho_0_i{iDim};
                    for i=1:3
                        if isnumeric(fun{i})
                            fun{i}=@(x)fun{i};
                        end
                    end
                end
                prjx=h(1)/2*base.projection(@(x)fun{1}((x+1)/2*h(1)+xyz0(1)),maxOrder,ngp);
                prjy=h(2)/2*base.projection(@(y)fun{2}((y+1)/2*h(2)+xyz0(2)),maxOrder,ngp);
                prjz=h(3)/2*base.projection(@(z)fun{3}((z+1)/2*h(3)+xyz0(3)),maxOrder,ngp);
                prj0=kron(kron(prjz,prjy),prjx);
                v0((iDim-1)*Nbasis+id)=v0((iDim-1)*Nbasis+id)+prj0;
            end
        end
    end
    
    % get u0 by solving linear system M*u0=v0 with v0=vec{(v_i,rho_0)}
    u0=zeros(Nbasis*dimRho,1);
    N_b=(1:Nbasis)';
    M=obj.M;
%     [massL,massU,massP,massQ,massR] = lu(M);
    for iDim=1:dimRho
%         u0(N_b+(iDim-1)*Nbasis)=massQ*(massU\(massL\(massP*(massR\v0(N_b+(iDim-1)*Nbasis)))));
        u0(N_b+(iDim-1)*Nbasis)=M\v0(N_b+(iDim-1)*Nbasis);
    end

end

function type=checkInitFormat(pars,rho_0)
    if iscell(rho_0) && isrow(rho_0) && length(rho_0)>1  % rho_0 is composed of many parts
        for i=1:length(rho_0)
            if ~pars.checkInitValFormat(rho_0{i})
                type=0; % format invalid
                return;  
            end
        end
        type=2;  % rho_0 is composed of many parts
        return;
    else
        if pars.checkInitValFormat(rho_0)
            type=1;  % rho_0 is composed of only one part
            return;  
        end
    end
    type=0; % format invalid
end
