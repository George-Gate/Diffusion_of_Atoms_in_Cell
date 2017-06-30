parSet3_new;  % get problemPars
% problemPars=PumpingDiffusionFEMSolver.ProblemPars();
% problemPars.T0=1;
FEMSolver=PumpingDiffusionFEMSolver.FEMsolver3D(problemPars);
FEMSolver.mesh.meshPars.Nz=8;
FEMSolver.mesh.meshPars.Nx=8;
FEMSolver.mesh.meshPars.Ny=8;
FEMSolver.simuPars.maxOrder=5;
FEMSolver.timeEvolution('parSet3, 8x8x8, maxOrder=5.mat');

%% plot result
Plotter=PumpingDiffusionFEMSolver.ResultVisualizer(FEMSolver);
Plotter.plotLines(1:3);
FEMSolver.plotMesh();
