parSet4_new;  % get problemPars
% problemPars=PumpingDiffusionFEMSolver.ProblemPars();
% problemPars.T0=1;
FEMSolver=PumpingDiffusionFEMSolver.FEMsolver3D(problemPars);
FEMSolver.mesh.meshPars.Nz=1;
FEMSolver.simuPars.maxOrder=6;
FEMSolver.timeEvolution();

%%
Plotter=PumpingDiffusionFEMSolver.ResultVisualizer(FEMSolver);
Plotter.plotLines(1);
FEMSolver.plotMesh();
