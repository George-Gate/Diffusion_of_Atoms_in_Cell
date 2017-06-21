function combineSolution( obj,sampleDim,lineID )
%Combine the solution and calc it on sampling points or let's say, mapping the solution to space-time domain.
%  [Usage]
%      sol = combineSolution( u, '3D', mesh, getNoByIxyz, K, Nbasis, dimRho, sampleRate )
%      sol = combineSolution( u, '1D', mesh, getNoByIxyz, K, Nbasis, dimRho, samplingLines )
%  [Input]
% --------------------- If sampleDim=='3D' ----------------------------------------------------------
%        u: the solution vector returned by ode45, each row is the solution for one time point.
%  sampleRate: the number of sampling points per surface per axis
%        K: the cut off of Lobatto Basis
%  mesh, getNoByIxyz: just pass those parameters
%
% --------------------- If sampleDim=='1D' ----------------------------------------------------------
%  samplingLines: the coordinate of sampling lines. Ns x 4 cell array, Ns is the number of sampling lines. 
%                 samplingLines{i,1:3}: xList/yList/zList of the sampling line, should be row or column vectors.
%                 samplingLines{i,4}:   the name or tag for this sampling line.
%
%  [Output]
% --------------------- If sampleDim=='3D' ----------------------------------------------------------
%    sol.rho : An (mesh.Ndomains x 1) cell array. sol.rho{i} is a (sampleRate x sampleRate x sampleRate x size(u,1) x dimRho) array containing 
%              the value of solution in domain i.
%    sol.xSample/.ySample/.zSample : The sampling grid in each domain. Is a mesh grid of [-1,1]^3.
%    sol.sampleDim = '3D'
%  
% --------------------- If sampleDim=='1D' ----------------------------------------------------------
%    sol.rho{i} : A (len_s x len_t x dimRho) array containing the value of solution.
%                 len_s is the length of sampling line, i.e. length(xList);
%                 len_t is the number of time points, i.e. size(u,1);
%    sol.xList{i}/.yList{i}/.zList{i} : Just store the input arguments, column vectors.
%    sol.tag{i} : the name/tag for current sampling line
%    sol.sampleDim = '1D'
%    
%

% lineID should be an integer

    switch sampleDim
        case '1D'
            if ~isempty(obj.sampleLineData{lineID})
                return;
            end
            % get parameters
            dimRho=obj.problemPars.dimRho;
            K=obj.simuPars.maxOrder;
            mesh=obj.mesh;
            base=obj.simuPars.baseFunHandle;
            Nbasis=base.Nbasis;
            getNoByIxyz=base.getNoByIxyz;
            % get current sample line
            samplingLines=obj.sampleLines{lineID};
            % transpose u, then each column of u stores the solution of one time point.
            u=transpose(obj.sol_u);
            % get coordinate list
            xList=samplingLines{1};  % these should be column vectors
            yList=samplingLines{2};
            zList=samplingLines{3};


       % ------------------------ combine solution --------------------------
            rho=zeros(length(xList),size(u,2),dimRho);

            for Did=1:mesh.Ndomains
                xyz=mesh.domains.xyz(:,:,Did);
                h=mesh.domains.h(:,Did);
                idRange= xyz(1,1)<=xList & xList<=xyz(1,2) & ...
                         xyz(2,1)<=yList & yList<=xyz(2,2) & ...
                         xyz(3,1)<=zList & zList<=xyz(3,2) ;
                %sol.rho(idRange,:,:)=0;
                len=sum(idRange);
                tmpSol=zeros(len,size(u,2),dimRho);
                % construct base value for this idRange
                coList=zeros(len,3);
                oneAxisBase=zeros(len,K,3);
                coList(:,1)=2*(xList(idRange)-xyz(1,1))/h(1)-1;
                coList(:,2)=2*(yList(idRange)-xyz(2,1))/h(2)-1;
                coList(:,3)=2*(zList(idRange)-xyz(3,1))/h(3)-1;
                for Aid=1:3
                    oneAxisBase(:,1:K,Aid)=base.funVal(coList(:,Aid),1:K); 
                end
                % enumerate basis  (This part should be optimized)
                for ix=1:K
                    for iy=1:K
                        for iz=1:K
                            iNo=getNoByIxyz(ix,iy,iz,Did);
                            if iNo==0; continue; end
                            % construct base value for this basis
                            baseValue=oneAxisBase(:,ix,1).*oneAxisBase(:,iy,2).*oneAxisBase(:,iz,3);
                            baseValue=repmat(baseValue,1,size(u,2));
                            % update tmpSol
                            for iDim=1:dimRho
                                tmpSol(:,:,iDim)=tmpSol(:,:,iDim)...
                                                +repmat(u(iNo+(iDim-1)*Nbasis,:),len,1).*baseValue;
                            end
                        end
                    end
                end
                % update sol
                rho(idRange,:,:)=tmpSol;
            end
            obj.sampleLineData{lineID}=rho;
            
        case '3D'
            error('This functionality is current not available.');
            sol.sampleDim=sampleDim;
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
            oneAxisBase=zeros(sampleRate,K);
            oneAxisBase(:,1)=(1-coList)/2;  % base @ x=-1
            oneAxisBase(:,2)=(1+coList)/2;  % base @ x=1
            oneAxisBase(:,3:K+2)=lobattoP_N(1:K-2,coList); % Lobatto base
            
            % --------------------- combine solution --------------------------
            % initial sol.rho
            sol.rho=cell(mesh.Ndomains,1);
            for Did=1:mesh.Ndomains
                sol.rho{Did}=zeros(sampleRate,sampleRate,sampleRate,size(u,2),dimRho);
            end
            
            % enumerate basis 
            for ix=1:K
                for iy=1:K
                    for iz=1:K
                        baseValueSet=0;
                        for Did=1:mesh.Ndomains
                            iNo=getNoByIxyz(ix,iy,iz,Did);
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

