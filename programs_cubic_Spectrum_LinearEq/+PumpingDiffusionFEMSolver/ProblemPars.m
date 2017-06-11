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
        rho_0;      % initial condition
        rho_b;      % boundary condition (first kind)
        drho_b;     % boundary condition (second kind)
        boundaryType; % boundary type, one of first, firstUniform, second, secondZero.
        anaSolForm; % analytical solution (if available)
    end
    
    properties(Dependent)
        dimRho;     % dimension of rho
        funP;       % 
        D;          % the diffusion coefficient D for the dimensionless PDE
        anaSol;     % analytical solution (if anaSolForm is set)
    end
    
    methods
        function obj=ProblemPars(filename)
            if nargin==0
                C=metaclass(obj);
                s=what(C.ContainingPackage.Name);
                filename=[s.path,'\DiffusionFEM.mat'];
            end
            diffPar=load(filename);
            obj.matC=diffPar.matC;
            obj.matD=diffPar.matD;
            obj.boundaryType='firstUniform';
            dimRho=size(obj.matD,1);
            obj.rho_b=ones(dimRho,1)/8/dimRho;
            obj.rho_0=ones(dimRho,1)/8/dimRho;
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
            obj.anaSolForm=fun;               % fun should be a function of obj....
        end
    end
    
end

