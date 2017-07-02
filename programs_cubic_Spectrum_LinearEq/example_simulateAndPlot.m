parSet3_new;  % get problemPars
% problemPars=PumpingDiffusionFEMSolver.ProblemPars();
% problemPars.T0=1;
FEMSolver=PumpingDiffusionFEMSolver.FEMsolver3D(problemPars);
FEMSolver.mesh.use_mex=0;
FEMSolver.mesh.meshPars.Nx=2;
FEMSolver.mesh.meshPars.Ny=2;
FEMSolver.mesh.meshPars.Nz=1;
FEMSolver.simuPars.maxOrder=12;
FEMSolver.timeEvolution('parSet3, 2x2x1, maxOrder=12.mat');

%% plot result
Plotter=PumpingDiffusionFEMSolver.ResultVisualizer(FEMSolver);
Plotter.plotLines(1);
Plotter.mesh.plotMesh();
