Nbasis=FEMSolver.baseFunction.Nbasis;
CB=FEMSolver.CB(1:Nbasis,1:Nbasis);

CB_correct=zeros(size(CB));
for iNo=1:size(CB,1)
    for kNo=1:size(CB,2)
        CB_correct(iNo,kNo)=calcCB(FEMSolver.baseFunction,FEMSolver.mesh,iNo,kNo);
    end
    disp(['iNo=',num2str(iNo),'/',num2str(Nbasis)]);
end

max(abs(CB(:)-CB_correct(:)))