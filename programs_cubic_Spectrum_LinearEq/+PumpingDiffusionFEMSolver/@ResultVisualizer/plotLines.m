function [  ] = plotLines(obj,lineIDs)
% Plot sampling lines
%   lineIDs - Sepcifies which sampling lines to plot

    num_lines=length(obj.sampleLines);
    if isempty(lineIDs) || ~isnumeric(lineIDs) || min(lineIDs(:))<1 || max(lineIDs(:))>num_lines
        error('Invalid Line No.');
    end
    lineIDs=reshape(lineIDs,1,numel(lineIDs));
    disp('Start plotting...... ');
    startT=tic;
    
    pars=obj.problemPars;
    maxOrder=obj.simuPars.maxOrder;
    dimRho=pars.dimRho;
    
% ----------- Plot Solution ---------------------------
    for Lid=lineIDs
        if isempty(obj.sampleLineData{Lid})
            % call combineSolution() method
            obj.combineSolution('1D',Lid);
        end
        xList=obj.sampleLines{Lid}{1};  % these should be column vectors
        yList=obj.sampleLines{Lid}{2};
        zList=obj.sampleLines{Lid}{3};
        tag=obj.sampleLines{Lid}{4};
        rho=obj.sampleLineData{Lid};
        % start plotting
        hf=figure('position',[0 0 1811 676],'name',['Sampling Line ',num2str(Lid),': ',tag]);
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
                plot(tList*pars.T0,8*sum(rho(:,:,iDim),3));
                ylabel(['$$L^3 \rho_',num2str(iDim),'$$'],'interpreter','latex');
            else  % 1D plot
                mesh(tList*pars.T0,xiList,8*sum(rho(:,:,iDim),3));
                view([0 90]);box on;
                colormap winter
                colorbar;
                ylabel('$$\xi$$','interpreter','latex');
                zlabel(['$$L^3 \rho_',num2str(iDim),'$$'],'interpreter','latex');
            end
            set(gca,'xlim',pars.T0*[min(tList(:)),max(tList(:))]);
            xlabel('time/s');
            title(['\rho_i, i=',num2str(iDim)]);
        end
        % the summation of all rho component
        subplot(3,4,dimRho+1);
        if isscalar(xList)
            % 0D plot
            plot(tList*pars.T0,(8*sum(rho(:,:,1:dimRho),3)-1));
            ylabel('$$L^3 \sum_{i}(\rho_i)-1$$','interpreter','latex');
        else  % 1D plot
            mesh(tList*pars.T0,xiList,(8*sum(rho(:,:,1:dimRho),3)-1));
            view([0 90]);box on;colormap winter;colorbar;
            ylabel('$$\xi$$','interpreter','latex');
            zlabel('$$L^3 \sum_{i}(\rho_i)-1$$','interpreter','latex');
        end
        set(gca,'xlim',pars.T0*[min(tList(:)),max(tList(:))]);
        xlabel('time/s');
        title('Normalization check.');
%         % average s_z
%         if isscalar(sol.xList{Lid})
%             % 0D plot
%             subplot(3,4,10);
%             len_s=size(sol.rho{Lid},1);len_t=size(sol.rho{Lid},2);
%             plot(tList*pars.T0,8*sum(sol.rho{Lid}(:,:,1:5).*repmat(reshape(-2:2,1,1,5),len_s,len_t,1),3));
%             box on;
%             set(gca,'xlim',pars.T0*[min(timeList),max(timeList)]);
%             xlabel('time/s');ylabel('Average F_z');
%             title('$$\langle F_z \rangle/\hbar \quad (\rho_1 - \rho_5)$$','interpreter','latex');
%             subplot(3,4,11);
%             plot(tList*pars.T0,8*sum(sol.rho{Lid}(:,:,6:8).*repmat(reshape(-1:1,1,1,3),len_s,len_t,1),3));
%             box on;
%             set(gca,'xlim',pars.T0*[min(timeList),max(timeList)]);
%             xlabel('time/s');ylabel('Average F_z');
%             title('$$\langle F_z \rangle/\hbar \quad (\rho_6 - \rho_8)$$','interpreter','latex');
%         else  % 1D plot
%             subplot(3,4,10);
%             len_s=size(sol.rho{Lid},1);len_t=size(sol.rho{Lid},2);
%             mesh(tList*pars.T0,xiList,8*sum(sol.rho{Lid}(:,:,1:5).*repmat(reshape(-2:2,1,1,5),len_s,len_t,1),3));
%             view([0 90]);box on;colormap winter;colorbar;
%             set(gca,'xlim',pars.T0*[min(timeList),max(timeList)]);
%             title('$$\langle F_z \rangle/\hbar \quad (\rho_1 - \rho_5)$$','interpreter','latex');
%             title('$$<F_z>/\hbar (\rho_1 -- \rho_5)$$','interpreter','latex');
%             subplot(3,4,11);
%             mesh(tList*pars.T0,xiList,8*sum(sol.rho{Lid}(:,:,6:8).*repmat(reshape(-1:1,1,1,3),len_s,len_t,1),3));
%             view([0 90]);box on;colormap winter;colorbar;
%             set(gca,'xlim',pars.T0*[min(timeList),max(timeList)]);
%             xlabel('time/s');ylabel('$$\xi$$','interpreter','latex');zlabel('Average F_z');
%             title('$$\langle F_z \rangle/\hbar \quad (\rho_6 - \rho_8)$$','interpreter','latex');
%         end
        % annotation
        annotation('textbox','String',...
           {['\rho(t),  K=',num2str(maxOrder),'  T_0=',num2str(pars.T0)],...
            ['Sampling Line ',num2str(Lid),': ',tag]}...
                  ,'Position',[0.75 0.21 0.1 0.1]);
        pause(0.1);
    end

    % show sample lines
    obj.showSampleLines(lineIDs);
    
    % display time consumption
    import PumpingDiffusionFEMSolver.library.sec2hms;
    disp(['Time used to plot sample lines: ',sec2hms(toc(startT))]);


end

