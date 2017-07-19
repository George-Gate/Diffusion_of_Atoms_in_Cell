noplot=1;
fileName='timeTest4';

NList=10000;
KpercentList=linspace(0.0005,0.01,2000);

testType='continuousAccess';
sparse_full_speedTest;
testType='jumpAccess';
sparse_full_speedTest;
% testType='randomAccess';
% sparse_full_speedTest;

