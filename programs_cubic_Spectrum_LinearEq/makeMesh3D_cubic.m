function  mesh  = makeMesh3D_cubic( xList, yList, zList )
%Generate a mesh for cubic
%   xList, yList, zList: vectors giving the partition of x, y, z axis

    % check if xList, yList and zList is ascending
    if min(diff(xList))<0 || min(diff(yList))<0 || min(diff(zList))<0
        warning('xList, yList and zList must be ascending, sort then automatically.');
        if min(diff(xList))<0
            xList=sort(xList);
        end
        if min(diff(yList))<0
            yList=sort(yList);
        end
        if min(diff(zList))<0
            zList=sort(zList);
        end
    end
        
    % calc element number
    Nx=length(xList);
    Ny=length(yList);
    Nz=length(zList);
    Nnodes=Nx*Ny*Nz;
    Nedges=(Nx-1)*Ny*Nz+Nx*(Ny-1)*Nz+Nx*Ny*(Nz-1);
    Nsurfaces=(Nx-1)*(Ny-1)*Nz+Nx*(Ny-1)*(Nz-1)+(Nx-1)*Ny*(Nz-1);
    Ndomains=(Nx-1)*(Ny-1)*(Nz-1);
    
    % def. mesh structure
    mesh.nodes.x=zeros(Nnodes,1);
    mesh.nodes.y=zeros(Nnodes,1);
    mesh.nodes.z=zeros(Nnodes,1);
    mesh.nodes.onB=zeros(Nnodes,1);         % .onB=1 if the node is on boundary; .onB=0 otherwise
    mesh.nodes.e=zeros(3,Nnodes);           % three edges whose .n(1,i) is current node, .e(1:3,i): edge parallel to x/y/z axis
    mesh.nodes.s=zeros(3,Nnodes);           % three surfaces whose .n(3,i) is current node, .s(1:3,i): surfaces that perpendicular to x/y/z axis
    mesh.edges.n=zeros(2,Nedges);           % .n(:,i): 2 endpoints, the x/y/z of n(1,i) is no larger than n(2,i)'s
    mesh.edges.onB=zeros(Nedges,1);
    mesh.surfaces.n=zeros(4,Nsurfaces);     % .n(:,i): 4 corner nodes, in quadrant 1~4. The third axis should be towards you (in right hand system). 
    mesh.surfaces.onB=zeros(Nsurfaces,1);
    mesh.domains.n=zeros(2,2,2,Ndomains);   % .n(dx,dy,dz,i): the node on relative position (dx,dy,dz); dx=1 for - , dx=2 for +
    mesh.domains.e=zeros(2,2,3,Ndomains);   % .e(dy,dz,1,i), .e(dz,dx,2,i), .e(dx,dy,3,i): the edge that parallel to x/y/z axis 
                                            % and on the relative position (dy,dz)/(dz,dx)/(dx,dy)
    mesh.domains.s=zeros(2,3,Ndomains);     % .s(dk,k,i): the surface that perpendicular to k axis and on the relative position dk       
    mesh.domains.h=zeros(3,Ndomains);       % .h(:,i): hx, hy and hz
    mesh.domains.xyz=zeros(3,2,Ndomains);   % .xyz(:,1,i): lower limit; .xyz(:,2,i): upper limit 
    mesh.Nnodes=Nnodes;
    mesh.Nedges=Nedges;
    mesh.Nsurfaces=Nsurfaces;
    mesh.Ndomains=Ndomains;    
    
    % set nodes, edges and surfaces
    topN=1;
    topEx=1;  topEy=topEx+(Nx-1)*Ny*Nz;  topEz=topEy+Nx*(Ny-1)*Nz;
    topSx=1;  topSy=topSx+Nx*(Ny-1)*(Nz-1);  topSz=topSy+(Nx-1)*Ny*(Nz-1);
    for nz=1:Nz
        for ny=1:Ny
            for nx=1:Nx
            % set nodes.x/.y/.z/.onB
                Nid=nx+(ny-1)*Nx+(nz-1)*Nx*Ny;
                if topN~=Nid
                    error('Nid not match!');
                end
                mesh.nodes.x(topN)=xList(nx);
                mesh.nodes.y(topN)=yList(ny);
                mesh.nodes.z(topN)=zList(nz);
                if (nx==1 || nx==Nx || ny==1 || ny==Ny || nz==1 || nz==Nz)
                    mesh.nodes.onB(topN)=1;
                end
            % set edges and nodes.e
                if (nx<Nx)  % edges parallel to x axis
                    mesh.edges.n(1:2,topEx)=[topN;topN+1];
                    mesh.nodes.e(1,topN)=topEx;
                    if (ny==1 || ny==Ny || nz==1 || nz==Nz)
                        mesh.edges.onB(topEx)=1;
                    end
                    topEx=topEx+1;
                end
                if (ny<Ny)   % edges parallel to y axis
                    mesh.edges.n(1:2,topEy)=[topN;topN+Nx];
                    mesh.nodes.e(2,topN)=topEy;
                    if (nx==1 || nx==Nx || nz==1 || nz==Nz)
                        mesh.edges.onB(topEy)=1;
                    end
                    topEy=topEy+1;
                end
                if (nz<Nz)   % edges parallel to z axis
                    mesh.edges.n(1:2,topEz)=[topN;topN+Nx*Ny];
                    mesh.nodes.e(3,topN)=topEz;
                    if (nx==1 || nx==Nx || ny==1 || ny==Ny)
                        mesh.edges.onB(topEz)=1;
                    end
                    topEz=topEz+1;
                end
            % set surfaces and nodes.s
                if (ny<Ny && nz<Nz)  % surface perp. to x axis
                    mesh.surfaces.n(1:4,topSx)=[topN+Nx+Nx*Ny;   topN+Nx*Ny;   topN;   topN+Nx];
                    mesh.nodes.s(1,topN)=topSx;
                    if (nx==1 || nx==Nx)
                        mesh.surfaces.onB(topSx)=1;
                    end
                    topSx=topSx+1;
                end
                if (nx<Nx && nz<Nz)  % surface perp. to y axis
                    mesh.surfaces.n(1:4,topSy)=[topN+1+Nx*Ny;   topN+1;   topN;   topN+Nx*Ny];
                    mesh.nodes.s(2,topN)=topSy;
                    if (ny==1 || ny==Ny)
                        mesh.surfaces.onB(topSy)=1;
                    end
                    topSy=topSy+1;
                end
                if (nx<Nx && ny<Ny)  % surface perp. to z axis
                    mesh.surfaces.n(1:4,topSz)=[topN+1+Nx;   topN+Nx;   topN;   topN+1];
                    mesh.nodes.s(3,topN)=topSz;
                    if (nz==1 || nz==Nz)
                        mesh.surfaces.onB(topSz)=1;
                    end
                    topSz=topSz+1;
                end
                % increase Nid to next one
                topN=topN+1;
            end
        end
    end
    
    % set domains
    topD=1;
    for nz=1:Nz-1
        for ny=1:Ny-1
            for nx=1:Nx-1
                Did=nx+(ny-1)*(Nx-1)+(nz-1)*(Nx-1)*(Ny-1);
                if topD~=Did
                    error('Did not match!');
                end
                Nid=nx+(ny-1)*Nx+(nz-1)*Nx*Ny;
                n=zeros(2,2,2);
                n(1:2,1:2,1)=[Nid,Nid+Nx;Nid+1,Nid+1+Nx];
                n(1:2,1:2,2)=[Nid,Nid+Nx;Nid+1,Nid+1+Nx]+Nx*Ny;
                mesh.domains.n(:,:,:,topD)=n;
                mesh.domains.e(:,:,1,topD)=reshape(mesh.nodes.e(1,[n(1,1,1),n(1,1,2);n(1,2,1),n(1,2,2)]),2,2);
                mesh.domains.e(:,:,2,topD)=reshape(mesh.nodes.e(2,[n(1,1,1),n(2,1,1);n(1,1,2),n(2,1,2)]),2,2);
                mesh.domains.e(:,:,3,topD)=reshape(mesh.nodes.e(3,[n(1,1,1),n(1,2,1);n(2,1,1),n(2,2,1)]),2,2);
                mesh.domains.s(:,1,topD)=reshape(mesh.nodes.s(1,n(:,1,1)),2,1);
                mesh.domains.s(:,2,topD)=reshape(mesh.nodes.s(2,n(1,:,1)),2,1);
                mesh.domains.s(:,3,topD)=reshape(mesh.nodes.s(3,n(1,1,:)),2,1);
                mesh.domains.xyz(1:3,1:2,topD)=[mesh.nodes.x(n(1,1,1)), mesh.nodes.x(n(2,2,2)) ;
                                                mesh.nodes.y(n(1,1,1)), mesh.nodes.y(n(2,2,2)) ;
                                                mesh.nodes.z(n(1,1,1)), mesh.nodes.z(n(2,2,2)) ];
                mesh.domains.h(:,topD)=mesh.domains.xyz(:,2,topD)-mesh.domains.xyz(:,1,topD);
                clear n
                topD=topD+1;
            end
        end
    end
    
end

