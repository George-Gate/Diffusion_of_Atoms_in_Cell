problem=PumpingDiffusionFEMSolver.ProblemPars;
% Assume that dimRho=8

% ------------ Example for vectors -------------------------

% constant for all boundary
problem.rho_b=ones(8,1)/8/8;   
% set boundary for each plane (x=-1, x=+1, y=-1, y=+1, z=-1, z=+1)
problem.drho_b={zeros(8,1);   % constant
                repmat({@(y,z)y.*z+1},8,1);  % 2D function
                repmat({{@(z)z,@(x)x+1}},8,1);  % separable 2D function
                zeros(8,1);
                {1;2;@(x,y)x.*y;4;5;{@(x)x,@(y)y+1};7;8}; % mixed
                zeros(8,1)};

% ------------ Example for matrices -------------------------

% constant matrix for all boundary
problem.robinA=ones(8);
% constant number for all boundary
problem.robinA=1;
% set boundary for each plane (x=-1, x=+1, y=-1, y=+1, z=-1, z=+1)
problem.robinA={{@(y,z)y+z};      % scalar 2D function
                {{@(y)2,@(z)1}};  % scalar separable 2D function
                 repmat({@(x,y)x+y},8);  % 2D function for each matrix elements
                 ones(8);  % constant matrix
                 1;        % constant scalar
                 {11,12,13,14,15,16,17,18;   % mixed form
                  21,22,23,24,25,26,27,28;
                  31,32,33,34,35,36,37,{@(x)x.^2;@(y)y};
                  41,42,43,44,45,46,47,@(x,y)x.*y+2;
                  51,52,53,54,55,56,57,58;
                  61,62,63,64,65,66,67,68;
                  71,72,73,74,75,76,77,78;
                  81,82,83,84,85,86,87,88;}};
% matrix-function separable form
problem.robinA={ones(8);          % constant matrix
                {@(y,z)y+z};      % scalar 2D function
                {ones(8),@(z,x)z+x};  % constant matrix times 2D function
                {ones(8),{@(z)z.*z,@(x)x.^2}};  % constant matrix times separable 2D function
                {{@(x)x,@(y)y}}; % scalar separable 2D function
                ones(8)};
           
% ------------ Example for initial condition -------------------------
problem.rho_0={1;    % constant
               2;    
               3;
               {4,5,6};  
               7;
               {@(x)x,@(y)y.^2,4};    % separable 3D function
               {@(x)x;@(y)y.^2;@(z)z}};
problem.rho_0=ones(8,1)/64;  % constant


