noplot=1;
fileName='timeTest3_CSRCpara';

NList=[60000:10000:140000];
KpercentList=linspace(0.0005,0.06,20);

testType='continuousAccess';
sparse_full_speedTest;
testType='jumpAccess';
sparse_full_speedTest;
testType='randomAccess';
sparse_full_speedTest;