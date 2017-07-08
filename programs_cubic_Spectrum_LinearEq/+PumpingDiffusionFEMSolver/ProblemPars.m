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
                    % It should be a function @(obj,t,x,y,z) with 
                    % obj of class @ProblemPars .
                    % For each (t,x,y,z) pair, anaSolForm should return
                    % a [dimRho x 1] column vector.
                    % anaSolForm can be designed to accept row vector input,
                    % i.e. given [1 x len] row vectors x,y,z and t as the input, 
                    % and return a [dimRho x len] matrix as output.
    end
    
    properties(Dependent)
        dimRho;     % dimension of rho
        funP;       % the spatial distribution of laser power
        D;          % the diffusion coefficient D for the dimensionless PDE
        anaSol;     % analytical solution (if anaSolForm is set), @(t,x,y,z)
    end
    
    methods
        function obj=ProblemPars(filename)
            if nargin==0
                C=metaclass(obj);
                s=what(C.ContainingPackage.Name);
                filename=fullfile(s.path,'DiffusionFEM.mat');
            end
            diffPar=load(filename);
            obj.matC=diffPar.matC;
            obj.matD=diffPar.matD;
            anaSolForm=@(obj,t,x,y,z)0;
            dimRho=obj.dimRho;
            obj.rho_0_ph=ones(dimRho,1)/dimRho/obj.L^3;
            % default boundary condition
            obj.boundaryType='second';
            obj.rho_b_ph=ones(dimRho,1)/dimRho/obj.L^3;
            obj.drho_b_ph=zeros(dimRho,1);
            obj.robinA_ph=0;
            obj.robinF_ph=zeros(dimRho,1);
        end
        
        function dimRho=get.dimRho(obj)
            if ~( size(obj.matD,1)==size(obj.matD,2) && size(obj.matC,1)==size(obj.matD,1) ...
               && size(obj.matC,1)==size(obj.matC,2))
                error('The size of matC and matD is invalid!');
            end
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
            anaSol=@(t,x,y,z)(obj.L/2)^3*obj.anaSolForm(obj,t*obj.T0,x*obj.L/2,y*obj.L/2,z*obj.L/2);
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
    end
% ============================= Initial state && Boundary condition ===========================================
        
    properties
        rho_0_ph;   % initial condition, should be a constant column vector 
                    % or a column cell vector that contains a separable 3D function
                    % in each element.
                    % New format: can be a ROW cell vector val with val{i} satisfying the above format.
                    % rho_0=val{1}+val{2}+val{3}+.....
        
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
        rho_b_ph;      % boundary condition (first type)
        drho_b_ph;     % boundary condition (second type)
        robinF_ph;     % boundary condition (robin type)
        robinA_ph;     % boundary condition (robin type)
        boundaryType;  % boundary type, one of first, second or robin.
    end

    properties(Dependent)
        rho_0;      % dimensionless form of rho_0_ph
        rho_b;      % dimensionless form of rho_b_ph
        drho_b;     % dimensionless form of drho_b_ph
        robinF;     % dimensionless form of robinF_ph
        robinA;     % dimensionless form of robinA_ph
    end

    methods
% ---------------------------- Set function for initial state -----------------------------------
        function obj=set.rho_0_ph(obj,val)
            errStr='The format of rho_0_ph is invalid.';
            if iscell(val) && isrow(val) && length(val)>1
                valid=true;
                for i=length(val)
                    if ~obj.checkInitValFormat(val{i})
                        valid=false;
                        break;
                    end
                end
            else
                valid=obj.checkInitValFormat(val);
            end
            if valid
                obj.rho_0_ph=val;
            else
                error(errStr);
            end
        end
        % Check whether a column vector val satisfies the initial value format.
        % Used by set.rho_0_ph() method and getInitialState() method of class @FEMSolver3D
        function valid=checkInitValFormat(obj,val)
            if isnumeric(val) && iscolumn(val) && length(val)==obj.dimRho
                valid=true; return;
            elseif iscell(val) && iscolumn(val) && length(val)==obj.dimRho
                % check form
                for i=1:obj.dimRho
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
                                valid=false; return;
                            end
                        end
                    else
                        valid=false; return;
                    end
                end
                valid=true; return;
            else
                valid=false; return;
            end
        end
% --------------------------------- Set functions of boundary conditions -------------------------------
        function obj=set.boundaryType(obj,val)
            validType={'first','second','robin'};
            if ismember(val,validType)
                obj.boundaryType=val;
            else
                error(['Unknow boundary type: ',val,'. Possible types are: ',strjoin(validType)]);
            end
        end
        function obj=set.rho_b_ph(obj,val)
            if obj.checkBoundaryForm_vector(val)
                obj.rho_b_ph=val;
            else
                error('The format of rho_b_ph is invalid.');
            end
        end
        function obj=set.drho_b_ph(obj,val)
            if obj.checkBoundaryForm_vector(val)
                obj.drho_b_ph=val;
            else
                error('The format of drho_b_ph is invalid.');
            end
        end
        function obj=set.robinF_ph(obj,val)
            if obj.checkBoundaryForm_vector(val)
                obj.robinF_ph=val;
            else
                error('The format of robinF_ph is invalid.');
            end
        end
        function obj=set.robinA_ph(obj,val)
            if obj.checkBoundaryForm_matrix(val)
                obj.robinA_ph=val;
            else
                error('The format of robinA_ph is invalid.');
            end
        end
        
% --------------------------------- Get functions of rho_0, rho_b, drho_b etc. -------------------------------
        function rho_0=get.rho_0(obj)
            % call setter to check format again in case that dimRho changed.
            obj.rho_0_ph=obj.rho_0_ph;  
            val=obj.rho_0_ph;
            if iscell(val) && isrow(val) && length(val)>1  % rho_0_ph is composed of many parts
                rho_0=cell(1,length(val));
                for i=1:length(val)
                    rho_0{i}=Nondimensionalization_init(val{i},obj.L,(obj.L/2)^3);
                end
            else
                rho_0=Nondimensionalization_init(val,obj.L,(obj.L/2)^3);
            end
        end
        function rho_b=get.rho_b(obj)
            % call setter to check format again in case that dimRho changed.
            obj.rho_b_ph=obj.rho_b_ph;  
            rho_b= Nondimensionalization_boundary(obj.rho_b_ph,  obj.L, (obj.L/2)^3, obj.dimRho);
        end
        function drho_b=get.drho_b(obj)
            % call setter to check format again in case that dimRho changed.
            obj.drho_b_ph=obj.drho_b_ph;  
            drho_b=Nondimensionalization_boundary(obj.drho_b_ph, obj.L, (obj.L/2)^4, obj.dimRho);
        end
        function robinA=get.robinA(obj)
            % call setter to check format again in case that dimRho changed.
            obj.robinA_ph=obj.robinA_ph;  
            robinA=Nondimensionalization_boundary(obj.robinA_ph, obj.L, (obj.L/2),   obj.dimRho);
        end
        function robinF=get.robinF(obj)
            % call setter to check format again in case that dimRho changed.
            obj.robinF_ph=obj.robinF_ph;  
            robinF=Nondimensionalization_boundary(obj.robinF_ph, obj.L, (obj.L/2)^4, obj.dimRho);
        end
    end
    
% -------------------------- Assistant Functions for format checking ----------------------------------------------
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

% ----------------- Assistant Functions for get functions of rho_0, rho_b, drho_b etc. ----------------------------------------------
% Get rho_0 from rho_0_ph
% Assume that rho_0_ph can pass the format check of initial value
function rho_0=Nondimensionalization_init(rho_0_ph,L,factor)
    dimRho=size(rho_0_ph,1);
    if isnumeric(rho_0_ph)       % constant init condition
        rho_0=rho_0_ph*factor;
    else                        % init condition if a function
        rho_0=cell(dimRho,1);
        for i=1:dimRho
            if isnumeric(rho_0_ph{i})
                rho_0{i}=rho_0_ph{i}*factor;
            else
                rho_0{i}=cell(1,3);
                factor3=factor^(1/3);
                for k=1:3
                    if isnumeric(rho_0_ph{i}{k})
                        rho_0{i}{k}=rho_0_ph{i}{k}*factor3;
                    else
                        rho_0{i}{k}=@(x)rho_0_ph{i}{k}(x*L/2)*factor3;
                    end
                end
            end
        end
    end
end

% Get rho_b from rho_b_ph, drho_b from drho_b_ph, robinA from robinA_ph, and robinF from robinF_ph
% Assume that bVal_ph can pass the format check of boundary condition
function bVal=Nondimensionalization_boundary(bVal_ph,L,factor,dimRho)
    if isnumeric(bVal_ph)    % bVal_ph is constant
        bVal=factor*bVal_ph;
    else                     % bVal_ph is a [6 x 1] cell vector
        bVal=cell(6,1);
        for iP=1:6
            if isnumeric(bVal_ph{iP})  % bVal_ph on this plane is constant
                bVal{iP}=factor*bVal_ph{iP};
            else             % bVal_ph on this plane is some spatial function, so it should be
                             % a cell array of size [1 x 1] or [dimRho x 1] or [dimRho x dimRho]
                             % or a special case for robinA: length 2 cell vector with bVal_ph0{1}
                             % a [dimRho x dimRho] numeric matrix and bVal_ph0{2} a 2D function
                bVal_ph0=bVal_ph{iP};
                bVal{iP}=cell(size(bVal_ph0));
                if isvector(bVal_ph0) && length(bVal_ph0)==2 ...  % special case for robinA
                    && size(bVal_ph0{1},1)==dimRho && size(bVal_ph0{1},2)==dimRho
                    bVal{iP}{1}=factor*bVal_ph0{1};
                    bVal{iP}{2}=convert2D(bVal_ph0{2},L,1);
                else                  % size=[1 x 1] or [dimRho x 1] or [dimRho x dimRho]
                    for ii=1:size(bVal_ph0,1)
                        for kk=1:size(bVal_ph0,2)
                            bVal{iP}{ii,kk}=convert2D(bVal_ph0{ii,kk},L,factor);
                        end
                    end
                end
            end
        end
    end
    
    % change the variables of a 2D function: f(x,y) -> factor*f(x*L/2,y*L/2)
    function fun2=convert2D(fun,L,factor)
        % fun should be either a numeric scalar or a 2D function handle or a length 2 cell vector
        % whose elements are either 1D function handles or numeric scalars
        if isnumeric(fun)
            fun2=factor*fun;
        elseif iscell(fun)
            fun2=cell(1,2);
            factor2=sqrt(factor);
            for i=1:2
                if isnumeric(fun{i})
                    fun2{i}=factor2*fun{i};
                else
                    fun2{i}=@(x)factor2*fun{i}(x*L/2);
                end
            end
        else
            fun2=@(x,y)factor*fun(x*L/2,y*L/2);
        end
    end
    
end