filename1='parSet1, 1x1x1, maxOrder=25.mat';
filename2='parSet1, 6x6x6, maxOrder=6.mat';
R1=load(filename1);
R2=load(filename2);


Plotter1=PumpingDiffusionFEMSolver.ResultVisualizer(R1.FEMResult);
Plotter2=PumpingDiffusionFEMSolver.ResultVisualizer(R2.FEMResult);

Plotter1.plotLines(5);
close;close;
Plotter2.plotLines(5);
close;close;

% change the SetAccess of sampleLineData to 'public' before running this script
Plotter1.sampleLineData{5}=Plotter1.sampleLineData{5}-Plotter2.sampleLineData{5};
Plotter1.plotLines(5);
