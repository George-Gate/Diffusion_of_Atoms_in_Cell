% Run time for K=12, nomesh
% Time to makeMesh: 0.06796 s
% Time to getCoeffs: 32.6723 s
% Time to Set Boundary: 0.031893 s
% Time to calc evolution: 3920.9115 s
% Time to combine solution: 145.779410 s.

% Run time for K=6, nomesh
% Time to makeMesh: 0.00998 s
% Time to getCoeffs: 0.84128 s
% Time to Set Boundary: 0.014168 s
% Time to calc evolution: 7.8871 s
% Time to save result to disk: 0.99902 s
% Time to combine solution: 16.9428 s


tic;
tid=[100];
combinedTime=reshape(tList(tid),length(tid),1);
sol=combineSolution( u(tid,:), '3D', mesh0, getNoByIxyz, K, Nbasis, dimRho, 100 );
toc;


%%
f=@(y)exp(-y.^2);
subid=4;
x=[-1,1];
xi=x(1);
h=x(2)-x(1);

I=h/2*integral(@(ksi)lobattoP_N(subid-1,ksi).*f( h/2*(ksi+1)+xi ) , -1,1)


% I=h/2*integral(@(ksi)lobattoP_N(subid1-1,ksi).*lobattoP_N(subid2-1,ksi).*f( h/2*(ksi+1)+xi ) , -1,1)
% (I-tmpIntGrid(subid1+1,subid2+1))/I

%%

 integral(@(x)funP{1}(x).*lobattoP_N(1,x),-1,1)...
*integral(@(y)funP{2}(y).*lobattoP_N(1,y),-1,1)...
*integral(@(z)funP{3}(z).*lobattoP_N(1,z),-1,1)...
*(matC(1,:)*rho_b)

 integral(@(x)lobattoP_N(1,x),-1,1)...
*integral(@(y)lobattoP_N(1,y),-1,1)...
*integral(@(z)lobattoP_N(1,z),-1,1)...
*(alpha*T0*matD(1,:)*rho_b)