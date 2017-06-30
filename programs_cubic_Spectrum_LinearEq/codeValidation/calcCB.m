function [ CB_ik ] = calcCB( base, mesh, iNo, kNo  )
%calcCB  Calc one element of matrix CB.
%   A very slow code, used for code validation
    
% ------------ Constant Definitions -----------------------------------------
    axis1=[2;3;1];      % axis1(dir)=the id for axis 1 if the plane direction is dir
    axis2=[3;1;2];      % axis2(dir)=the id for axis 2 if the plane direction is dir
    axis3=[1;2;3];      % axis3=dir, the plane is orthonal to axis 3
    xaxis=[3;2;1];      % xaxis(dir)=the axis number of x axis when plane direction is dir
    yaxis=[1;3;2];      
    zaxis=[2;1;3];


    ngp=base.maxOrder*2+100;
    getNoByIxyz=base.getNoByIxyz;
    
    [ix,iy,iz,iDid]=ind2sub(size(getNoByIxyz),find(getNoByIxyz==iNo));
    [kx,ky,kz,kDid]=ind2sub(size(getNoByIxyz),find(getNoByIxyz==kNo));
    CB_ik=0;
    sub_i=[ix,iy,iz];
    sub_k=[kx,ky,kz];
    for i=1:length(iDid)
        Did=iDid(i);
        if ~ismember(Did,kDid)
            continue;
            % each boundary surface only belongs to one domain.
            % So in order to have non-zero value on the same surface,
            % base iNo and kNo should have non-zero value in the domain
            % that the boundary surface belongs to simultaneously. 
        end
        k=find(kDid==Did);
        % the base now is iNo -> (ix(i),iy(i),iz(i),Did);  kNo -> (kx(k),ky(k),kz(k),Did)
        for dk=1:2
            for dir=1:3
                % skip non-boundary syrfaces
                Sid=mesh.domains.s(dk,dir,Did);
                if ~mesh.surfaces.onB(Sid)
                    continue;
                end
                % find integral area by surfaces.n[1] and surfaces.n[3]
                n1=mesh.surfaces.n(1,Sid);
                n3=mesh.surfaces.n(3,Sid);
                xyz=[mesh.nodes.x([n1;n3]),mesh.nodes.y([n1;n3]),mesh.nodes.z([n1;n3])];
                range1=[xyz(2,axis1(dir));xyz(1,axis1(dir))];
                range2=[xyz(2,axis2(dir));xyz(1,axis2(dir))];
                planePos=xyz(1,axis3(dir));
                if xyz(1,axis3(dir))~=xyz(2,axis3(dir))
                    error('Coordinate check failed.');
                end
                h1=range1(2)-range1(1);h2=range2(2)-range2(1);
                
                int1=h1/2*base.projection(@(x)base.funVal( x,sub_i(i,axis1(dir)) ), sub_k(k,axis1(dir)), ngp);
                int1=int1(end);
                int2=h2/2*base.projection(@(x)base.funVal( x,sub_i(i,axis2(dir)) ), sub_k(k,axis2(dir)), ngp);
                int2=int2(end);
                co3=base.funVal(planePos, sub_i(i,axis3(dir)) )*base.funFirstDerivative(planePos, sub_k(k,axis3(dir)) ) ...
                   +base.funVal(planePos, sub_k(k,axis3(dir)) )*base.funFirstDerivative(planePos, sub_i(i,axis3(dir)) );
                co3=co3*sign(planePos);
                CB_ik=CB_ik+int1*int2*co3;
            end
        end
    end


end

