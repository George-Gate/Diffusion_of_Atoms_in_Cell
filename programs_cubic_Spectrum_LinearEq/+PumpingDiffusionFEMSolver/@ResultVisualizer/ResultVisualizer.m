classdef ResultVisualizer < handle
    %ResultVisualizer  Provide visualization for the result of FEM solver
    
    properties(SetAccess=private)
        sampleLines;       % sample line list
        sampleLineData;    % data cache for sample line plotting
        
        % information of FEM solution
        sol_u;sol_t;
        solSampleRate;   % The sampling rate set in @FEMsolver3D object for time evolution
        mesh;
        simuPars;
        problemPars;
    end
    
    
    methods(Access=private)
        % Combine the solution and calc it on sampling points
        combineSolution(obj,sampleDim,lineID);
    end
    
    methods
        % Plot the selected sample lines
        plotLines(obj,lineIDs);
    end
    
    
    methods
        function obj=ResultVisualizer(FEMResult)
            % initialize sampleLines & sampleLineData
            obj.sampleLines=defaultLines();
            obj.sampleLineData=repmat({{}},length(obj.sampleLines),1);
            
            if exist('FEMResult','var')
                obj.setFEMResult(FEMResult);
            end
        end
        
        % Add a new sample line to sampleline list, the sample line should be inside [-1,1]^3
        function addSampleLine(obj,xList,yList,zList,label)
            if ischar(label) && isnumeric(xList) && isnumeric(yList) && isnumeric(zList) ...
                    && isvector(xList) && isvector(yList) && isvector(zList) ...
                    && length(xList)==length(yList) && length(yList)==length(zList)
                
                num_lines=length(obj.sampleLines);
                len=length(xList);
                % reshape xList, tList and zList to column vector
                xList=reshape(xList,len,1);yList=reshape(yList,len,1);zList=reshape(zList,len,1);
                obj.sampleLines{num_lines+1}={xList,yList,zList,label};
                obj.sampleLineData{num_lines+1}={};
            else
                error('Input format invalid. xList, yList and zList should be three numerical vectors with the same length. label should be a char array.');
            end
        end
        
        % Print out the sample line list
        function listSampleLines(obj)
            num_lines=length(obj.sampleLines);
            if num_lines==0
                disp('No sample line.');
            else
                fprintf('%-8s|      %-20s|      %-25s\n','Line No.','(x,y,z)','label');
                disp('--------+--------------------------+-------------------------------');
                for i=1:num_lines
                    if length(obj.sampleLines{i}{1})==1
                        xyz=[obj.sampleLines{i}{1:3}];
                        coStr=['(',num2str(xyz(1),2),', ',num2str(xyz(2),2),', ',num2str(xyz(3),2),')'];
                    else
                        coStr=[num2str(length(obj.sampleLines{i}{1})),' pts.'];
                    end
                    fprintf('   %-5i| %-25s| %-30s\n',i,coStr,obj.sampleLines{i}{4});
                end
            end
        end
        
        % Show specific sample lines in a figure
        function showSampleLines(obj,lineIDs)
            if nargin<2
                lineIDs=1:length(obj.sampleLines);
            end
            if isempty(lineIDs)
                disp('No sample line.');
            else
                if ~isnumeric(lineIDs) || min(lineIDs(:))<1 || max(lineIDs(:))>length(obj.sampleLines)
                    error('Invalid Line No.');
                end
                lineIDs=reshape(lineIDs,1,numel(lineIDs));
            % -------------------- Plot Sampling Lines ---------------------------------
                figure();
                tag={};
                for Lid=lineIDs
                    xList=obj.sampleLines{Lid}{1};  % these should be column vectors
                    yList=obj.sampleLines{Lid}{2};
                    zList=obj.sampleLines{Lid}{3};
                    tag=[tag,{obj.sampleLines{Lid}{4}}]; %#ok<AGROW>
                    if isscalar(xList)
                        scatter3(xList,yList,zList);
                    else
                        plot3(xList,yList,zList);
                    end
                    hold on;
                end
                legend(tag{:});
                title('Sampling Lines');box on;
                xlabel('x');ylabel('y');zlabel('z');
                set(gca,'xlim',[-1,1],'ylim',[-1,1],'zlim',[-1,1]);
            end
        end
        
        % Delete some sample lines from sample line list
        function delSampleLines(obj,lineIDs)
            num_lines=length(obj.sampleLines);
            if ~isnumeric(lineIDs) || min(lineIDs(:))<1 || max(lineIDs(:))>num_lines
                error('Invalid Line No.');
            end
            remainList=[];
            for i=1:num_lines
                if ~ismember(i,lineIDs)
                    remainList=[remainList;i]; %#ok<AGROW>
                end
            end
            obj.sampleLines={obj.sampleLines{remainList}}';
            obj.sampleLineData={obj.sampleLineData{remainList}}';
        end
        
        % Load the result of FEM solver
        function setFEMResult(obj,FEMResult)
            if isa(FEMResult,'PumpingDiffusionFEMSolver.FEMsolver3D')
                FEMResult=FEMResult.FEMResult;
            end
            if isstruct(FEMResult) && isfield(FEMResult,'sol_u') && isfield(FEMResult,'sol_t') ...
                    && isfield(FEMResult,'meshPars') && isfield(FEMResult,'simuPars') && isfield(FEMResult,'problemPars') ...
                    && isfield(FEMResult,'sampleRate')
                obj.sol_u=FEMResult.sol_u;
                obj.sol_t=FEMResult.sol_t;
                obj.mesh=PumpingDiffusionFEMSolver.Mesh3D(FEMResult.meshPars);
                obj.solSampleRate=FEMResult.sampleRate;
                obj.simuPars=FEMResult.simuPars;
                obj.simuPars.baseFunHandle=obj.simuPars.baseFunHandle.copy();   % make a deep copy to avoid affecting the FEMsolver3D
                obj.problemPars=FEMResult.problemPars;
                % re-generate mesh and base numbering
                obj.mesh.makeMesh();
                if ~ (obj.simuPars.baseFunHandle.meshPars==obj.mesh.meshPars_current ...
                        && obj.simuPars.baseFunHandle.maxOrder==obj.simuPars.maxOrder)
                    obj.simuPars.baseFunHandle.numbering(obj.mesh,obj.simuPars.maxOrder);
                end
                % clean cache data
                obj.sampleLineData=repmat({{}},length(obj.sampleLines),1);
            else
                error(['The type of input is not supported. Please give an @FEMsolver3D object or the FEMResult structure generated by this object.', ...
                       ' If the input is an @FEMsolver3D object, please ensure that the timeEvolution() method is called before.']);
            end
        end
        
        
% ------------------ Setters & Getters -------------------------------------------------------------------------------
        function set.simuPars(obj,val)
            if isa(val,'PumpingDiffusionFEMSolver.SimuPars')
                obj.simuPars=val;
            else
                error('simuPars should be an object of class @PumpingDiffusionFEMSolver.SimuPars');
            end
        end
        function set.problemPars(obj,val)
            if isa(val,'PumpingDiffusionFEMSolver.ProblemPars')
                obj.problemPars=val;
            else
                error('simuPars should be an object of class @PumpingDiffusionFEMSolver.ProblemPars');
            end
        end
        function set.mesh(obj,val)
            if isa(val,'PumpingDiffusionFEMSolver.Mesh3D')
                obj.mesh=val;
            else
                error('mesh should be an object of class @PumpingDiffusionFEMSolver.Mesh3D');
            end
        end
        
        
    end
    
end


% return default samplelines as an example
function lines=defaultLines()
    lines={{linspace(-1,1,200),zeros(1,200),zeros(1,200),'xAxis'};
           {zeros(1,200),linspace(-1,1,200),zeros(1,200),'yAxis'};
           {zeros(1,200),zeros(1,200),linspace(-1,1,200),'zAxis'};
           {linspace(-1,1,200),zeros(1,200),-0.99*ones(1,200),'(y=0,z=-0.99)'};
           {linspace(-1,1,200),linspace(-1,1,200),linspace(-1,1,200),'(-1,-1,-1)->(1,1,1)'};
           {0,0,0,'centerPoint'};
           {0.5,0.5,0.5,'Point:(0.5,0.5,0.5)'};
           {1,1,1,'Point:(1,1,1)'};};
end
