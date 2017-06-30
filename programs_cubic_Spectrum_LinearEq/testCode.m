CMG=PumpingDiffusionFEMSolver.CoeffMatrixGenerator();
CMG.simuPars.maxOrder=6;
CMG.mesh.meshPars.zMajor=[-1;0;0.5;1];
CMG.mesh.meshPars.Nz=[1,1,1];
CMG.mesh.meshPars.yMajor=[-1;0.1;1];
CMG.mesh.meshPars.Ny=[1,1];
CMG.problemPars.rho_0={1;
                       2;
                       3;
                       {4,5,6};
                       7;
                       {@(x)1,@(y)3,4};
                       {@(x)x;3;@(z)z}};
% CMG.problemPars.boundaryType='second';
% CMG.problemPars.drho_b={ones(8,1);
%                        {1;2;3;4;5;@(x,y)x+y;7;8};
%                        {1;2;3;4;5;{@(x)x,@(z)z+1};7;8};
%                        ones(8,1);
%                        repmat({@(x,y)x+y},8,1);
%                        ones(8,1)};
CMG.problemPars.boundaryType='robin';
CMG.problemPars.robinF={ones(8,1);
                        {1;2;3;4;5;@(x,y)x+y;7;8};
                        {1;2;3;4;5;{@(x)x,@(z)z+1};7;8};
                        ones(8,1);
                        repmat({@(x,y)x+y},8,1);
                        ones(8,1)};
CMG.problemPars.robinA={ones(8);          % constant matrix
                        {@(y,z)y+z};      % scalar 2D function
                        {ones(8),@(z,x)z+x};  % constant matrix times 2D function
                        {ones(8),{@(z)z.*z,@(x)x.^2}};  % constant matrix times separable 2D function
                        {{@(x)x,@(y)y}}; % scalar separable 2D function
                        ones(8)};
CMG.genCoeffs

%%
robinA={ones(8);          % constant matrix
                        {@(y,z)y+z};      % scalar 2D function
                        {ones(8),@(z,x)z+x};  % constant matrix times 2D function
                        {ones(8),{@(z)z.*z,@(x)x.^2}};  % constant matrix times separable 2D function
                        {{@(x)x,@(y)y}}; % scalar separable 2D function
                        ones(8)};
robinB=robinA;