% testType='randomAccess';
% testType='continuousAccess';
% testType='jumpAccess';

% fileName='timeTest1';

% NList=[1000];
% KpercentList=[0.001:0.001:0.009,0.01:0.01:0.09,0.1:0.1:0.5];

% timeOfSparse(iK,iN,1): average; timeOfSparse(iK,iN,2): std dev
timeOfSparse=zeros(length(KpercentList),length(NList),2);
timeOfFull=timeOfSparse;


for iN=1:length(NList)
    N=NList(iN);
    iList0=(1:N)';
    jList0=(1:N)';
    vList0=rand(N,1);
    tNstart=tic;
    for iK=1:length(KpercentList)
        K=max(1,round(N*KpercentList(iK)));

        vList=repmat(vList0,K,1);
        iList=repmat(iList0,K,1);
        if strcmp(testType,'jumpAccess')
            jList=repmat(jList0,K,1);
        else
            jList=repelem(jList0,K,1);
        end
        summer_s=0; summer_f=0; avgN=20;
        summer_s2=0; summer_f2=0;
        for i=1:avgN
            if strcmp(testType,'randomAccess')
                id=randperm(K*N);
                jList=jList(id);
                clear id;
            end

            tic;
            S=sparse(iList,jList,vList,N,N);
            t=toc;   clear S;
            summer_s=summer_s+t;  summer_s2=summer_s2+t*t;
            tic;
            F=vec2full(iList,jList,vList,N,N);
            t=toc;  clear F;
            summer_f=summer_f+t;  summer_f2=summer_f2+t*t;
        end
        timeOfSparse(iK,iN,1)=summer_s/avgN;
        timeOfFull(iK,iN,1)=summer_f/avgN;
        timeOfSparse(iK,iN,2)=sqrt(summer_s2/(avgN-1)-(summer_s^2/avgN/(avgN-1)));
        timeOfFull(iK,iN,2)=sqrt(summer_f2/(avgN-1)-(summer_f^2/avgN/(avgN-1)));
        save([fileName,'_',testType,'.mat'],'timeOfSparse','timeOfFull','NList','KpercentList','testType','avgN');
    end
    disp(['Test for N=',num2str(N),' Finished. Save result. Time:',num2str(toc(tNstart),3),'s']);
end


%% plot
if ~exist('noplot','var') || ~noplot
    figure();
    KpercentList=reshape(KpercentList,length(KpercentList),1);
    NList=reshape(NList,length(NList),1);
    subplot(3,1,1);
    for iN=1:length(NList)
        relErrS=timeOfSparse(:,iN,2)./timeOfSparse(:,iN,1);
        relErrF=timeOfFull(:,iN,2)./timeOfFull(:,iN,1);
        absErr=timeOfSparse(:,iN,1)./timeOfFull(:,iN,1).*sqrt(relErrS.^2+relErrF.^2);
%         errorbar(KpercentList,timeOfSparse(:,iN,1)./timeOfFull(:,iN,1),absErr);hold on
        plot(KpercentList,timeOfSparse(:,iN,1)./timeOfFull(:,iN,1));hold on
    end
    ylabel('Relative time consumption. t_{Sparse}/t_{Full}');
    xlabel('Density of matrix');
    title(['testType=',testType,'  avgN=',num2str(avgN)]);
    legend(num2str(NList));

    subplot(3,1,2);
    for iN=1:length(NList)
        errorbar(KpercentList,timeOfSparse(:,iN,1),timeOfSparse(:,iN,2));hold on
    end
    ylabel('Time consumption of Sparse()');
    xlabel('Density of matrix');
    legend(num2str(NList));

    subplot(3,1,3);
    for iN=1:length(NList)
        errorbar(KpercentList,timeOfFull(:,iN,1),timeOfFull(:,iN,2));hold on
    end
    ylabel('Time consumption of Full()');
    xlabel('Density of matrix');
    legend(num2str(NList));
end
