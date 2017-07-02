function status=timeEvolutionOutputFunction(obj,t,y,flag)
%timeEvolutionOutputFunction  The Output Function for ode solvers
% Functionality: Output evolution progress and save result during evolution
    status=0;
    switch flag
        case 'done'
            % do cleaning
            obj.lastOutputTime=uint64(0);
            obj.lastReportedT=0;
            obj.lastSaveTime=uint64(0);
            obj.lastSavedT=0;
            obj.tyRecord=[];
            % output message
            disp(['[',datestr(datetime,'mmm-dd HH:MM:SS'),'] Time evolution completed!']);
        case []
            % record result
            top=obj.tyRecord.top;
            lenT=length(t);
            obj.tyRecord.t(top:top+lenT-1)=t;
            obj.tyRecord.y(:,top:top+lenT-1)=y;
            obj.tyRecord.top=top+lenT;
            % disp simulation progress
            if t(end)-obj.lastReportedT>0.1 || toc(obj.lastOutputTime) > 300
                estRemainTime=toc(obj.lastOutputTime)/(t(end)-obj.lastReportedT)*(1-t(end));
                import PumpingDiffusionFEMSolver.library.sec2hms;
                disp(['[',datestr(datetime,'mmm-dd HH:MM:SS'),num2str(100*t(end),'] %4.1f'),num2str(t(end),'%%  t=%4.3f'),'    Remaining: ',sec2hms(estRemainTime)]);
                obj.lastOutputTime=tic;
                obj.lastReportedT=t(end);
            end
            % save intermediate results to file
            if (toc(obj.lastSaveTime)>3600)   % save every 1 hour
                id=find(obj.tyRecord.t>obj.lastSavedT);
                tyRecord.t=obj.tyRecord.t(id);
                tyRecord.y=obj.tyRecord.y(:,id); %#ok<STRNU>
                filename=['tyRecord_t0=',num2str(obj.lastSavedT),'.mat'];
                save(filename,'tyRecord');
                disp(['[',datestr(datetime,'mmm-dd HH:MM:SS'),'] Calculated solution saved to ''',filename,'''.']);
                obj.lastSaveTime=tic;
                obj.lastSavedT=t;
            end
        case 'init'
            % init flags
            obj.lastOutputTime=tic;
            obj.lastReportedT=0;
            obj.lastSaveTime=tic;
            obj.lastSavedT=0;
            % init obj.tyRecord
            recordLen=obj.sampleRate*obj.problemPars.T0+1;
            obj.tyRecord.t=zeros(1,recordLen);
            obj.tyRecord.y=zeros(obj.problemPars.dimRho*obj.baseFunction.Nbasis,recordLen);
            obj.tyRecord.top=1;
            % save parameters
            filename=obj.saveResultToFile('parameters.mat');
            disp(['[',datestr(datetime,'mmm-dd HH:MM:SS'),'] Parameters saved to ''',filename,'''.']);
            % output message
            disp(['[',datestr(datetime,'mmm-dd HH:MM:SS'),'] First time step!']);
    end
end