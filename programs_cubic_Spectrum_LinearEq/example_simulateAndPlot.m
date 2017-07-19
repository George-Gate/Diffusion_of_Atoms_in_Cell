parSet1;
% parSet5_new;  % get problemPars
% problemPars=PumpingDiffusionFEMSolver.ProblemPars();
% problemPars.T0=1;
filename='parSet1, 1x1x1, maxOrder=10.mat';
disp(filename);
FEMSolver=PumpingDiffusionFEMSolver.FEMsolver3D(problemPars);
FEMSolver.mesh.use_mex=0;
FEMSolver.mesh.meshPars.Nx=1;
FEMSolver.mesh.meshPars.Ny=1;
FEMSolver.mesh.meshPars.Nz=1;
FEMSolver.simuPars.maxOrder=10;
% FEMSolver.sampleRate=1000;
FEMSolver.timeEvolution(filename);
FEMResult=FEMSolver.FEMResult;

%% plot result
Plotter=PumpingDiffusionFEMSolver.ResultVisualizer(FEMResult);
Plotter.plotLines(5);
% Plotter.compareWithAnaSol_line(5);
% Plotter.mesh.plotMesh();
