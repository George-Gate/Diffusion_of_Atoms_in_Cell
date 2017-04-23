function sol = combineSolution( u, sampleDim, mesh, getNoByIxyz_a, K, Nbasis, dimRho, varargin )
%Combine the solution and calc it on sampling points
%  [Usage]
%      sol = combineSolution( u, '3D', mesh, getNoByIxyz_a, K, Nbasis, dimRho, sampleRate )
%  [Input]
% --------------------- If sampleDim=='3D' ----------------------------------------------------------
%        u: the solution vector returned by ode45, each row is the solution for one time point.
%  sampleRate: the number of sampling points per surface per axis
%        K: the cut off of Lobatto Basis
%  mesh, getNoByIxyz_a: just pass those parameters
%
%  [Output]
%    sol: An (Nsurfaces x 1) cell array. Each element is a (sampleRate x sampleRate) matrix containing 
%         the value of solution in each surface.
    switch sampleDim
        case '3D'
            % parse arguments
            sampleRate=varargin{1};
            % transpose u, then each column of u stores the solution of one time point.
            u=transpose(u);  
            % get sampling grid
            [xSample,ySample,zSample]=meshgrid(linspace(-1,1,sampleRate));
            sol.xSample=xSample;
            sol.ySample=ySample;
            sol.zSample=zSample;
            % calc one axis base value
            coList=reshape(xSample(1,:,1),sampleRate,1);
            oneAxisBase=zeros(sampleRate,K+2);
            oneAxisBase(:,1)=(1-coList)/2;  % base @ x=-1
            oneAxisBase(:,2)=(1+coList)/2;  % base @ x=1
            oneAxisBase(:,3:K+2)=lobattoP_N(1:K,coList); % Lobatto base
            
            % --------------------- combine solution --------------------------
            % initial sol.rho
            sol.rho=cell(mesh.Ndomains,1);
            for Did=1:mesh.Ndomains
                sol.rho{Did}=zeros(sampleRate,sampleRate,sampleRate,size(u,2),dimRho);
            end

            % enumerate basis 
            for ix=0:K+1
                for iy=0:K+1
                    for iz=0:K+1
                        baseValueSet=0;
                        for Did=1:mesh.Ndomains
                            iNo=getNoByIxyz_a(ix+1,iy+1,iz+1,Did);
                            if iNo==0; continue; end
                            % construct base value for this basis
                            if ~baseValueSet
                                xPart=reshape(oneAxisBase(:,ix+1),1,sampleRate,1);
                                yPart=reshape(oneAxisBase(:,iy+1),sampleRate,1,1);
                                zPart=reshape(oneAxisBase(:,iz+1),1,1,sampleRate);
                                xPart=repmat(xPart,sampleRate,1,sampleRate);
                                yPart=repmat(yPart,1,sampleRate,sampleRate);
                                zPart=repmat(zPart,sampleRate,sampleRate,1);
                                baseValue=xPart.*yPart.*zPart;
                                baseValueSet=1;
                            end
                            for it=1:size(u,2)
                                for iDim=1:dimRho
                                    sol.rho{Did}(:,:,:,it,iDim)=sol.rho{Did}(:,:,:,it,iDim)...
                                                                +u(iNo+(iDim-1)*Nbasis,it)*baseValue;
                                end
                            end
                        end
                    end
                end
            end

        otherwise
            error(['Unknow sampleDim: ',sampleDim]);
    end
end

