function [  ] = compareWithNumResult_line(obj1,lineIDs1,obj2,lineIDs2)
% Compare the solution of obj1 and obj2 on specific sampling lines
%   obj1, obj2 - @ResultVisualizer objects
%   lineIDs1, lineIDs2 - Sepcifies which sampling lines to compare

    % check lineIDs
    num_lines1=length(obj1.sampleLines);
    if ~isnumeric(lineIDs1) || min(lineIDs1(:))<1 || max(lineIDs1(:))>num_lines1
        error('Invalid Line No.');
    end
    lineIDs1=reshape(lineIDs1,1,numel(lineIDs1));
    
    num_lines2=length(obj2.sampleLines);
    if ~isnumeric(lineIDs2) || min(lineIDs2(:))<1 || max(lineIDs2(:))>num_lines2
        error('Invalid Line No.');
    end
    lineIDs2=reshape(lineIDs2,1,numel(lineIDs2));
    
    % check pars compatibility
    if ~(obj1.problemPars==obj2.problemPars)
        error('Cannot compare these two results because their problemPars are different.');
    end
    
    disp('Start plotting...... ');
    startT=tic;
    
    pars=obj1.problemPars;
 %%!!!   maxOrder=obj1.simuPars.maxOrder;
    dimRho=pars.dimRho;
    
% ----------- Plot numSol - anaSol ---------------------------
    for Lid=lineIDs1
        if isempty(obj1.sampleLineData{Lid})
            % call combineSolution() method
            obj1.combineSolution('1D',Lid);
        end
        xList=obj1.sampleLines{Lid}{1};  % these should be column vectors
        yList=obj1.sampleLines{Lid}{2};
        zList=obj1.sampleLines{Lid}{3};
        tag=obj1.sampleLines{Lid}{4};
        rho_num=obj1.sampleLineData{Lid};
        % calc rho_ana
        rho_ana=calcAnaSol(pars,obj1.sol_t,xList,yList,zList,dimRho);
        
        % start plotting
        hf=figure('position',[0 0 1811 676],'name',['numSol-anaSol. Sampling Line ',num2str(Lid),': ',tag]);
        if isscalar(xList)
            % 0D plot
            tList=reshape(obj1.sol_t,length(obj1.sol_t),1);  % tList should be column vector
        else  % 1D plot
            [tList,xiList]=meshgrid(obj1.sol_t,linspace(0,1,length(xList)));
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

