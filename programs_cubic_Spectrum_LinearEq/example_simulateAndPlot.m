anaSol_2;
% parSet5_new;  % get problemPars
% problemPars=PumpingDiffusionFEMSolver.ProblemPars();
% problemPars.T0=1;
filename='anaSol_2, 6x6x6, maxOrder=6_3.mat';
disp(filename);
FEMSolver=PumpingDiffusionFEMSolver.FEMsolver3D(problemPars);
FEMSolver.mesh.use_mex=0;
FEMSolver.mesh.meshPars.Nx=6;
FEMSolver.mesh.meshPars.Ny=6;
FEMSolver.mesh.meshPars.Nz=6;
FEMSolver.simuPars.maxOrder=6;
FEMSolver.timeEvolution(filename);
FEMResult=FEMSolver.FEMResult;

%% plot result
% Plotter=PumpingDiffusionFEMSolver.ResultVisualizer(FEMResult);
% % Plotter.plotLines(2);
% Plotter.compareWithAnaSol_line(1:7);
% % Plotter.mesh.plotMesh();
