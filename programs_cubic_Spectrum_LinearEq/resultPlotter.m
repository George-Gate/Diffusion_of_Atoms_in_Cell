function [  ] = resultPlotter( sol, pars, timeList )
% Plot solutions returned be combineSolution()
%   sol - returned by combineSolution()
%   pars - parameters used to solve equation
%   samplingLines - samplingLines 

fprintf('Plot...... ');tic;
% plot1D
if strcmp(sol.sampleDim,'1D')
    for Lid=1:length(sol.tag)
        hf=figure('position',[1426 0 1811 676],'name',['Sampling Line ',num2str(Lid),': ',sol.tag{Lid}]);
        if isscalar(sol.xList{Lid})
            % 0D plot
            tList=timeList;
        else  % 1D plot
            [tList,xiList]=meshgrid(timeList,linspace(0,1,length(sol.xList{Lid})));
        end
        for rhoComponent=1:8
            subplot(3,4,rhoComponent);
            if isscalar(sol.xList{Lid})
                % 0D plot
                plot(tList*pars.T0,sum(sol.rho{Lid}(:,:,rhoComponent),3));
                ylabel('Normalized \rho-\rho_0');
            else  % 1D plot
                mesh(tList*pars.T0,xiList,sum(sol.rho{Lid}(:,:,rhoComponent),3));
                view([0 90]);box on;
                colormap winter
                colorbar;
                ylabel('$$\xi$$','interpreter','latex');zlabel('Normalized \rho-\rho_0');
            end
            set(gca,'xlim',pars.T0*[min(timeList),max(timeList)]);
            xlabel('time/s');
            title(['\rho_i, i=',num2str(rhoComponent)]);
        end
        % the summation of all rho component
        subplot(3,4,9);
        if isscalar(sol.xList{Lid})
            % 0D plot
            plot(tList*pars.T0,sum(sol.rho{Lid}(:,:,1:8),3));
            ylabel('Normalized \rho-\rho_0');
        else  % 1D plot
            mesh(tList*pars.T0,xiList,sum(sol.rho{Lid}(:,:,1:8),3));
            view([0 90]);box on;colormap winter;colorbar;
            ylabel('$$\xi$$','interpreter','latex');zlabel('Normalized \rho-\rho_0');
        end
        set(gca,'xlim',pars.T0*[min(timeList),max(timeList)]);
        xlabel('time/s');
        title('sum(\rho_i)');
        % average s_z
%         subplot(3,4,10);
%         len_s=size(sol.rho{Lid},1);len_t=size(sol.rho{Lid},2);
%         mesh(tList*pars.T0,xiList,sum(sol.rho{Lid}(:,:,1:8).*repmat(reshape(-3.5:3.5,1,1,8),len_s,len_t,1),3));
%         view([0 90]);box on;colormap winter;colorbar;
%         set(gca,'xlim',pars.T0*[min(timeList),max(timeList)]);
%         xlabel('time/s');ylabel('$$\xi$$','interpreter','latex');zlabel('Average S_z');
%         title('<S_z>');
        % annotation
        annotation('textbox','String',...
           {['\rho(t)-\rho(0),  K=',num2str(pars.K),'  T_0=',num2str(pars.T0)],...
            ['Sampling Line ',num2str(Lid),': ',sol.tag{Lid}]}...
                  ,'Position',[0.55 0.21 0.1 0.1]);
        pause(0.01);
    end
    figure();
    for Lid=1:length(sol.tag)
        if isscalar(sol.xList{Lid})
            scatter3(sol.xList{Lid},sol.yList{Lid},sol.zList{Lid});
        else
            plot3(sol.xList{Lid},sol.yList{Lid},sol.zList{Lid});
        end
        hold on;
    end
    legend(sol.tag{:});
    title('Sampling Lines');box on;
    xlabel('x');ylabel('y');zlabel('z');
% plot3D
elseif strcmp(sol.sampleDim,'3D')
    hf=figure('position',[1426 82 1811 533]);
    tid=1;
    for rhoComponent=1:8
        subplot(2,4,rhoComponent);
        hold on;
        for Did=1:mesh0.Ndomains
            midPos=sum(mesh0.domains.xyz(:,:,Did),2)/2;
            domW=mesh0.domains.h(:,Did)/2;
            xGrid=sol.xSample*domW(1)+midPos(1);
            yGrid=sol.ySample*domW(2)+midPos(2);
            zGrid=sol.zSample*domW(3)+midPos(3);
            hs=slice(xGrid,yGrid,zGrid,sum(sol.rho{Did}(:,:,:,tid,rhoComponent),5),midPos(1),midPos(2),midPos(3));
            for i=1:length(hs)
                hs(i).EdgeColor='none';
            end
        end
        view([26 6]);
        colormap winter
        colorbar;
        title(['\rho_',num2str(rhoComponent)]);
        xlabel('x');ylabel('y');zlabel('z');
    end
    titleTextBox=uicontrol(hf,'style','text');
    annotation('textbox','String',...
        {['\rho(t)-\rho(0),  K=',num2str(K)],...
         ['T=',num2str(pars.T0*timeList(tid)),' s']}...
              ,'Position',[0.02 0.8 0.1 0.1]);
end
disp([num2str(toc),' s']);


end

