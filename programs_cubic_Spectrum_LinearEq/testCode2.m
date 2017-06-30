invMH=FEMSolver.MM\FEMSolver.H;
invMF=FEMSolver.MM\FEMSolver.F;
eigList=eig(invMH);
[min(real(eigList)) max(real(eigList))]