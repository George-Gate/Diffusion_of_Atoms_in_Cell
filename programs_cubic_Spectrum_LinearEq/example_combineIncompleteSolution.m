%% load and combine result
% tyRecord file list, in time increasing order
fileList={'tyRecord_t0=0.mat';
          'tyRecord_t0=0.048.mat';
          'tyRecord_t0=0.188.mat'};
load parameters.mat
FEMResult.sol_t=[];
FEMResult.sol_u=[];
for fID=1:length(fileList)
    load(fileList{fID});
    FEMResult.sol_t=[FEMResult.sol_t,tyRecord.t];
    FEMResult.sol_u=[FEMResult.sol_u,tyRecord.y];
end
FEMResult.sol_u=FEMResult.sol_u.';


%% plot result
Plotter=PumpingDiffusionFEMSolver.ResultVisualizer(FEMResult);
Plotter.plotLines(1:3);
Plotter.mesh.plotMesh();
