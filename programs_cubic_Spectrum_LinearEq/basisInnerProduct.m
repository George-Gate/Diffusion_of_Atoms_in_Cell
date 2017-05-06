function I = basisInnerProduct( fun, getNoByIxyz, K, Nbasis, mesh, ngp)
%Calc the inner product of all base given function 
%
%       I = basisInnerProduct( fun, getNoByIxyz, K, Nbasis, mesh, ngp)
%
%   For Bilinear+Lobatto basis
%   If fun is a constant function (a number), the computation will be very fast.
%   The returned value I is a ( Nbasis x 1 ) vector. I(:) is fun*(v_i,1), i=1,2,3,...,Nbasis
%
%   If fun is a (3 x 1) cell array and fun{1/2/3} is a function_handle with one input argument or scalar number.
%   This program will assume fun(x,y,z)=fun{1}(x)*fun{2}(y)*fun{3}(z).
%   The returned value I is a column vector of length Nbasis, giving the integrals (v_i, fun(x,y,z) ) 
%
%
%
%
    if isnumeric(fun) && isscalar(fun)
        I=zeros(Nbasis,1);
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
                        I(No)=I(No)+Iv*fun;
                    end
                end
            end
        end
    elseif iscell(fun) && isvector(fun) && length(fun)==3 
    
        % --------------- Generate gauss points and weights --------------
        % gp_x: gauss points; gw: weights
        [ gp_x, gw ] = getGaussPts( ngp );
        % ----------------- Evaluate basis ----------------------------
        phi_num=zeros(ngp,K+2);
        phi_num(:,1)=(1-gp_x)/2;
        phi_num(:,2)=(1+gp_x)/2;
        phi_num(:,3:K+2)=lobattoP_N(1:K,gp_x);
        % --------------- Calc Integrals ---------------------------------
        I=ones(Nbasis,1);   % I(No)=(v_No,fun(x,y,z))
        for Did=1:mesh.Ndomains
            for dir=1:3
                if isa(fun{dir},'function_handle')
                    intRange=mesh.domains.xyz(dir,:,Did)';
                    P_num=gw.*fun{dir}( (gp_x+1)/2*(intRange(2)-intRange(1))+intRange(1) );
                    tmpIntList=( intRange(2)-intRange(1) )/2*phi_num'*P_num; % tmpIntList(l+1)=(phi_l,funP);
                elseif isnumeric(fun{dir}) && isscalar(fun{dir})
                    tmpIntList=fun{dir}*[1;1;-2/3;zeros(K-1,1)];
                else
                    error(['Unsupported type of fun{',num2str(dir),'}']);
                end
                for l=0:K+1   % times tmpIntList into I
                    if dir==1
                        NoList=getNoByIxyz(l+1,:,:,Did);
                    elseif dir==2
                        NoList=getNoByIxyz(:,l+1,:,Did);
                    else
                        NoList=getNoByIxyz(:,:,l+1,Did);
                    end
                    NoList=NoList(NoList>0);
                    I(NoList)=I(NoList)*tmpIntList(l+1);
                end
            end
        end
    else
        error('Unsupported type of fun');
    end
end

