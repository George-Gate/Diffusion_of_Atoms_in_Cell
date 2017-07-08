function [  ] = compareWithAnaSol_line(obj,lineIDs)
% Plot sampling lines
%   lineIDs - Sepcifies which sampling lines to plot

    num_lines=length(obj.sampleLines);
    if ~isnumeric(lineIDs) || min(lineIDs(:))<1 || max(lineIDs(:))>num_lines
        error('Invalid Line No.');
    end
    lineIDs=reshape(lineIDs,1,numel(lineIDs));
    disp('Start plotting...... ');
    startT=tic;
    
    pars=obj.problemPars;
    maxOrder=obj.simuPars.maxOrder;
    dimRho=pars.dimRho;
    
% ----------- Plot numSol - anaSol ---------------------------
    for Lid=lineIDs
        if isempty(obj.sampleLineData{Lid})
            % call combineSolution() method
            obj.combineSolution('1D',Lid);
        end
        xList=obj.sampleLines{Lid}{1};  % these should be column vectors
        yList=obj.sampleLines{Lid}{2};
        zList=obj.sampleLines{Lid}{3};
        tag=obj.sampleLines{Lid}{4};
        rho_num=obj.sampleLineData{Lid};
        % calc rho_ana
        rho_ana=calcAnaSol(pars,obj.sol_t,xList,yList,zList,dimRho);
        
        % start plotting
        hf=figure('position',[0 0 1811 676],'name',['numSol-anaSol. Sampling Line ',num2str(Lid),': ',tag]);
        if isscalar(xList)
            % 0D plot
            tList=reshape(obj.sol_t,length(obj.sol_t),1);  % tList should be column vector
        else  % 1D plot
            [tList,xiList]=meshgrid(obj.sol_t,linspace(0,1,length(xList)));
        end
        for iDim=1:dimRho
            subplot(3,4,iDim);
            if isscalar(xList)
                % 0D plot
                plot(tList*pars.T0,8*sum(rho_num(:,:,iDim)-rho_ana(:,:,iDim),3));
                ylabel(['$$L^3 \left( \rho_{\mathrm{N},',num2str(iDim),'} - \rho_{\mathrm{A},',num2str(iDim),'} \right) $$'],'interpreter','latex');
            else  % 1D plot
                mesh(tList*pars.T0,xiList,8*sum((rho_num(:,:,iDim)-rho_ana(:,:,iDim)),3));
                view([0 90]);box on;
                colormap winter
                colorbar;
                ylabel('$$\xi$$','interpreter','latex');
                zlabel(['$$L^3 \left( \rho_{\mathrm{N},',num2str(iDim),'} - \rho_{\mathrm{A},',num2str(iDim),'} \right) $$'],'interpreter','latex');
            end
            set(gca,'xlim',pars.T0*[min(tList(:)),max(tList(:))]);
            xlabel('time/s');
            title(['\rho_i, i=',num2str(iDim)]);
        end

        % annotation
        annotation('textbox','String',...
           {['\rho(t),  K=',num2str(maxOrder),'  T_0=',num2str(pars.T0)],...
            ['Sampling Line ',num2str(Lid),': ',tag]}...
                  ,'Position',[0.75 0.21 0.1 0.1]);
        pause(0.1);
    end
    
    
    % display time consumption
    import PumpingDiffusionFEMSolver.library.sec2hms;
    disp(['Time used to plot sample lines: ',sec2hms(toc(startT))]);


end


function rho_ana=calcAnaSol(pars,t,x,y,z,dimRho)
    if ~(isvector(x) && isvector(y) && isvector(z) && isvector(t) ...
         && length(x)==length(y) && length(y)==length(z))
        error('Dimension check of x, y, z and t failed.');
    end
    len=length(x);lenT=length(t);
    x=reshape(x,1,len);y=reshape(y,1,len);z=reshape(z,1,len);
    t=reshape(t,1,lenT);
   
    anaFun=pars.anaSol;
    % check if anaFun can accept vector input
    vecInput=true;
    try
        anaFun([0,0.5],[0,0.5],[0,0.5],[0,0.5]);
        anaFun([0,0.5,1],[0,0.5,1],[0,0.5,1],[0,0.5,1]);
    catch
        vecInput=false;
    end
    rho_ana=zeros(len,lenT,dimRho);
    if ~vecInput
        for i=1:len
            for j=1:lenT
                tmp=anaFun(t(j),x(i),y(i),z(i));
                rho_ana(i,j,:)=reshape(tmp,1,1,dimRho);
            end
        end
    else
        if len>lenT    % there are less time points than spatial points
                       % so use x, y, z as vector inputs for each time point
            for j=1:lenT
                tmp=anaFun(t(j),x,y,z);
                rho_ana(:,j,:)=permute(tmp,[2 3 1]);
            end
        else
            for i=1:len
                tmp=anaFun(t,x(i),y(i),z(i));
                rho_ana(i,:,:)=permute(tmp,[3 2 1]);
            end
        end
    end
    
end