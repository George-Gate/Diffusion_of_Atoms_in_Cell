classdef ProblemPars
    %ProblemPars  Provide PDE parameters for @FEMsolver3D
    
    properties
        matC;       % pumping efficiency
        matD;       % relaxation rate
        alpha=0.5;  % 0.5/s
        w=0.2;      % 0.2cm
        P0=1;       % 1/s
        D_ph=0.1;   % 0.1 cm^2/s
        T0=0.5;     % evolution time, in unit of seconds
        L=1;        % length of cubic cell, in unit of cm
        anaSolForm; % analytical solution (if available)
        rho_0;      % initial condition, should be a constant column vector 
                    % or a column cell vector that contains three 1D functions
                    % in each element.
        
        %-------------[Format of boundary condition]------------------
        % 1. A [dimRho x 1] constant vector, saying that the boundary value if independent of space 
        % 2. A [6 x 1] cell array, each element can be either a [dimRho x 1] constant vector or
        %    a [dimRho x 1] cell array that contain 2-D function handles which returns a scalar.
        %    The 2-D function handle should be something like 
        %                      @(x,y)... or @(y,z)... or @(z,x)... 
        %                   or {@(x)...;@(y)...} or {@(y)...;@(z)...} or {@(z)...;@(x)...}
        %    The [6 x 1] cell array provides the boundary value for x=-1, x=1, y=-1, y=1, z=-1, z=1 plane.
        %##  For the case of the matrix A in robin boundary condition, the size [dimRho x 1] above is 
        %    replaced by [dimRho x dimRho] or [1 x 1].
        %### For robinA, there is a new matrix-function separable form.
        %    Each element of the [6 x 1] cell array can be a length 2 cell vector like the following:
        %            {rand(dimRho),@(z,x)z+x}
        %            {rand(dimRho),{@(z)z.*z,@(x)x.^2}}
        % !!! The boudary condition should be time independent.
        rho_b;      % boundary condition (first type)
        drho_b;     % boundary condition (second type)
        robinF;     % boundary condition (robin type)
        robinA;     % boundary condition (robin type)
        boundaryType; % boundary type, one of first, second, robin.
    end
    
    properties(Dependent)
        dimRho;     % dimension of rho
        funP;       % the spatial distribution of laser power
        D;          % the diffusion coefficient D for the dimensionless PDE
        anaSol;     % analytical solution (if anaSolForm is set)
    end
    
    methods
        function obj=ProblemPars(filename)
            if nargin==0
                C=metaclass(obj);
                s=what(C.ContainingPackage.Name);
                if ispc
                    filename=[s.path,'\DiffusionFEM.mat'];
                elseif isunix
                    filename=[s.path,'/DiffusionFEM.mat'];
                end
            end
            diffPar=load(filename);
            obj.matC=diffPar.matC;
            obj.matD=diffPar.matD;
            anaSolForm=@(obj,t,x,y,z)0;
            dimRho=obj.dimRho;
            obj.rho_0=ones(dimRho,1)/8/dimRho;
            % default boundary condition
            obj.boundaryType='first';
            obj.rho_b=ones(dimRho,1)/8/dimRho;
            obj.drho_b=zeros(dimRho,1);
            obj.robinA=0;
            obj.robinF=zeros(dimRho,1);
        end
        
        function dimRho=get.dimRho(obj)
            dimRho=size(obj.matD,1);
        end
        function D=get.D(obj)
            D=4*obj.T0/obj.L/obj.L*obj.D_ph;
        end
        function funP=get.funP(obj)
            funP={@(x)exp( -x.*x * (obj.L/2/obj.w)^2 );...
                  @(y)exp( -y.*y * (obj.L/2/obj.w)^2 );...
                  @(z)obj.T0*obj.P0};
        end
        function anaSol=get.anaSol(obj) % !!!!! this part not finished
            anaSol=@(t,x,y,z)obj.anaSolForm(obj,t,x,y,z);
        end
        function obj=set.anaSolForm(obj,fun)  % !!!!! this part not finished
            if isa(fun,'function_handle') && nargin(fun)==5
                obj.anaSolForm=fun;               % fun should be a function of obj....
            else
                error('fun should be a function handle: fun=@(obj,t,x,y,z)...');
            end
        end
% ====================== Definition of eq operator ============================
        function yes=eq(o1,o2)
            yes=true;
            if ~isa(o1,'PumpingDiffusionFEMSolver.ProblemPars') || ...
                    ~isa(o2,'PumpingDiffusionFEMSolver.ProblemPars')
                yes=false;
                return;
            end
            if ~( isequal(o1.matC,o2.matC) && isequal(o1.matD,o2.matD) ...
               && o1.alpha==o2.alpha && o1.w==o2.w ...
               && o1.P0==o2.P0 && o1.D_ph==o2.D_ph ...
               && o1.T0==o2.T0 && o1.L==o2.L ...
               && isequal(o1.rho_0,o2.rho_0) && strcmp(o1.boundaryType,o2.boundaryType))
                yes=false;
                return;
            end
            switch o1.boundaryType
                case 'first'
                    if ~( isequal(o1.rho_b,o2.rho_b) )
                        yes=false;
                        return;
                    end
                case 'second'
                    if ~( isequal(o1.drho_b,o2.drho_b) )
                        yes=false;
                        return;
                    end
                case 'robin'
                    if ~( isequal(o1.robinF,o2.robinF) && isequal(o1.robinA,o2.robinA) )
                        yes=false;
                        return;
                    end
            end
        end
        

% ====================== set function for initial state ============================
        function obj=set.rho_0(obj,val)
            if isnumeric(val) && iscolumn(val)
                obj.rho_0=val;
            elseif iscell(val) && iscolumn(val)
                for i=1:length(val)
                    if isnumeric(val{i}) && isscalar(val{i})
                        % good
                    elseif iscell(val{i}) && isvector(val{i}) && length(val{i})==3
                        % check whether there are three 1D function
                        tmp=val{i};
                        for j=1:3
                            if isnumeric(tmp{j})
                                % good
                            elseif isa(tmp{j},'function_handle') && nargin(tmp{j})==1
                                % good
                            else
                                error('The format of rho_0 is invalid.');
                            end
                        end
                    else
                        error('The format of rho_0 is invalid.');
                    end
                end
            else
                error('The format of rho_0 is invalid.');
            end
        end
% ====================== set functions of boundary conditions ============================
        function obj=set.boundaryType(obj,val)
            validType={'first','second','robin'};
            if ismember(val,validType)
                obj.boundaryType=val;
            else
                error(['Unknow boundary type: ',val,'. Possible types are: ',strjoin(validType)]);
            end
        end
        function obj=set.rho_b(obj,val)
            if obj.checkBoundaryForm_vector(val)
                obj.rho_b=val;
            else
                error('The format of rho_b is invalid.');
            end
        end
        function obj=set.drho_b(obj,val)
            if obj.checkBoundaryForm_vector(val)
                obj.drho_b=val;
            else
                error('The format of drho_b is invalid.');
            end
        end
        function obj=set.robinF(obj,val)
            if obj.checkBoundaryForm_vector(val)
                obj.robinF=val;
            else
                error('The format of robinF is invalid.');
            end
        end
        function obj=set.robinA(obj,val)
            if obj.checkBoundaryForm_matrix(val)
                obj.robinA=val;
            else
                error('The format of robinA is invalid.');
            end
        end
% ======================================================================================  
    end
    
    methods(Hidden,Access=private)
        % This function defines the format of boundary condition
        % obj.rho_b, obj.drho_b, obj.robinF
        function valid=checkBoundaryForm_vector(obj,val)
            if isnumeric(val)                           % format 1
                valid=isNumericColVector(val,obj.dimRho);
            elseif iscell(val) && isvector(val) && length(val)==6    % format 2
                valid=true;
                for i=1:6              % check the form of each component
                    tmp=val{i};
                    if iscell(tmp)
                        tmp=cellfun2num(tmp);
                    end
                    if ~isNumericColVector(tmp,obj.dimRho)
                        valid=false;
                        break;
                    end
                end
            else
                valid=false;
            end
        end
        % this function defines the format of boundary condition 
        % obj.robinA
        function valid=checkBoundaryForm_matrix(obj,val)
            if isnumeric(val)                           % format 1
                valid=isNumericMatrix(val,obj.dimRho);
            elseif iscell(val) && isvector(val) && length(val)==6    % format 2
                valid=true;
                for i=1:6              % check the form of each cell component
                    tmp=val{i};
                    
                    % matrix-function separable form
                    if iscell(tmp) && isvector(tmp) && (length(tmp)==2 || length(tmp)==3) ...
                           && isNumericMatrix(tmp{1},obj.dimRho)
                        
                        fun_num=cellfun2num({tmp{2:end}});
                        if isscalar(fun_num)
                            tmp=tmp{1}*fun_num;
                        else
                            tmp=[];
                        end
                    % element-wise form
                    elseif iscell(tmp)
                        tmp=cellfun2num(tmp);
                    end
                    if ~isNumericMatrix(tmp,obj.dimRho)
                        valid=false;
                        break;
                    end
                end
            else
                valid=false;
            end
        end
    end
    
end

function num=cellfun2num(funArray)
    num=zeros(size(funArray));
    for i=1:numel(funArray)
        fun=funArray{i};
        if isnumeric(fun) && isscalar(fun)
            num(i)=fun;
        elseif isa(fun,'function_handle') && nargin(fun)==2
            num(i)=fun(0,0);
        elseif iscell(fun) && isvector(fun) && isa(fun{1},'function_handle') && nargin(fun{1})==1 ...
                                            && isa(fun{2},'function_handle') && nargin(fun{2})==1
            num(i)=fun{1}(0)*fun{2}(0);
        else
            num=[];
            break;
        end
    end
end

function yes=isNumericColVector(val,len)
    if isnumeric(val) && iscolumn(val) && length(val)==len
        yes=true;
    else
        yes=false;
    end
end

function yes=isNumericMatrix(val,len)
    if isnumeric(val) && ismatrix(val) && size(val,1)==len && size(val,2)==len
        yes=true;
    elseif isnumeric(val) && isscalar(val)
        yes=true;
    else
        yes=false;
    end
end
