%% load and combine result
% tyRecord file list, in time increasing order
% folder path
folder=fullfile('Simulation_Results','parSet3, 2x2x1, maxOrder=14');
% get file list
fileObj=dir(fullfile(folder,'tyRecord_t0=*.mat'));
fileList={fileObj.name}';
% load data from files
load parameters.mat
FEMResult.sol_t=0;
FEMResult.sol_u=FEMResult.u0;
for fID=1:length(fileList)
    load(fullfile(folder,fileList{fID}));
    FEMResult.sol_t=[FEMResult.sol_t,tyRecord.t];
    FEMResult.sol_u=[FEMResult.sol_u,tyRecord.y];
end
% sort by time
[FEMResult.sol_t,I]=sort(FEMResult.sol_t);
FEMResult.sol_u=FEMResult.sol_u(:,I);
% transpose
FEMResult.sol_u=FEMResult.sol_u.';


%% plot result
Plotter=PumpingDiffusionFEMSolver.ResultVisualizer(FEMResult);
Plotter.plotLines(1);
Plotter.mesh.plotMesh();
