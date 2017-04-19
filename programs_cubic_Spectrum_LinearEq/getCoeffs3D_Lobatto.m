function [ Mvec, Svec, MPintegralList, QintegralList  ] = getCoeff_Lobatto( boundaryType, mesh, K )
% Generate the coeff matrices for a box area.
%    Do not support non constant first kind boundary condition so far. 
%    boundaryType: 'first', 'second', 'firstConstant', 'secondZero'
%    if boundaryType == 'firstConstant', we will give a matrix that solves zero boundary problem. i.e. our solution is (u-c) where c is the boundary value.
%
%    K: the cut off of Lobatto polynomials
%
%   MPintegralList: a list of integral to be calculated. matrix MP depends on these integrals.
%     MPintegralList(i,k).ndim
%                        .ub/.lb
%                        .fun :     @(weightFun,x,y,z,...)
%      We suppose there is a function that can calc the integral: 
%             int(@(x,y,z...).fun(weightFun,x,y,z,...),lb,ub);
%      Where weightFun is a ndim-D function given by user
%   
%   QintegralList: similar to MPintegralList, its value give vecQ
%
%
%
%


% ==================== read mesh ================================
if strcmp(boundaryType,'firstConstant') || strcmp(boundaryType,'first')
    boundaryBasis=0;
elseif strcmp(boundaryType,'secondZero') || strcmp(boundaryType,'second')
    boundaryBasis=1;
else
    error(['Unknow boundary type: ',boundaryType]);
end
% set fun2No and No2fun
if (boundaryBasis)
    NBnode=mesh.Nnodes;
    NBedge=mesh.Nedges*K;
    NBface=mesh.Nsurfaces*K*K;
else
    NBnode=mesh.Nnodes-sum(mesh.nodes.onB);
    NBedge=(mesh.Nedges-sum(mesh.edges.onB))*K;
    NBface=(mesh.Nsurfaces-sum(mesh.surfaces.onB))*K*K;
end
NBdomain=mesh.Ndomains*K*K*K;
Nbasis=NBnode+NBedge+NBface+NBdomain;

% array def.
fun2No.node=zeros(mesh.Nnodes,1);
fun2No.edge=zeros(K+1,mesh.Nedges);              % subid=1: linear basis;  subid=2~K+1: Lobatto basis
fun2No.face=zeros(K+1,K+1,mesh.Nsurfaces);       % fun2No.face(subid_x,subid_y,Fid)/(subid_y,subid_z,Fid)/(subid_z,subid_x,Fid) for different faceDir
fun2No.domain=zeros(K+1,K+1,K+1,mesh.Ndomains);  % fun2No.domain(subid_x, subid_y, subid_z, Did)
No2fun.name=char(zeros(Nbasis,1));   % possible values: 'n' for node; 'e' for edge; 'f' for face; 'd' for domain
No2fun.objid=zeros(Nbasis,1);
No2fun.subid=zeros(Nbasis,3);        % subid=1: linear basis;  subid=2~K+1: Lobatto basis; subid=[subid_x;subid_y;subid_z]
edgeDir=zeros(mesh.Nedges,1);        % 1/2/3: parallel to x/y/z axis
faceDir=zeros(mesh.Nsurfaces,1);     % 1/2/3: perpendicular to x/y/z axis


% node mode
No2fun.name(1:NBnode)='n';
No2fun.subid(1:NBnode,1:3)=1;
if boundaryBasis
    No2fun.objid(1:NBnode)=(1:mesh.Nnodes)';
else
    No2fun.objid(1:NBnode)=find(mesh.nodes.onB==0);
end
fun2No.nodal(No2fun.objid(1:NBnode))=(1:Nnodal)';

% edge mode
No2fun.name(NBnode+1:NBnode+NBedge)='e';
topNo=NBnode+1;
for Nid=1:mesh.Nnodes
    for iEid=1:3
        Eid=mesh.nodes.e(iEid,Nid);
        edgeDir(Eid)=iEid;
        if (mesh.edges.onB(Eid) && ~boundaryBasis)
            % skip boundary elements if don't need boundary basis
            continue;
        end
        for subid=2:K+1
            fun2No.edge(subid,Eid)=topNo;
            No2fun.objid(topNo)=Eid;
            No2fun.subid(topNo,:)=1;
            No2fun.subid(topNo,iEid)=subid;
            topNo=topNo+1;
        end
    end
end
if ~( No2fun.name(topNo)=='' &&  No2fun.name(topNo-1)=='e' )
    error ('Basis number check error!');
end

% face mode
No2fun.name(topNo:topNo+NBface)='f';
for Nid=1:mesh.Nnodes
    for iSid=1:3
        Sid=mesh.nodes.s(iSid,Nid);
        faceDir(Sid)=iSid;
        if (mesh.surfaces.onB(Sid) && ~boundaryBasis)
            % skip boundary elements if don't need boundary basis
            continue;
        end
        for subid1=2:K+1
            for subid2=2:K+1
                fun2No.face(subid1,subid2,Sid)=topNo;
                No2fun.objid(topNo)=Sid;
                if iSid==1      % face perp. to x axis   subid1 -> subid_y; subid2 -> subid_z
                    No2fun.subid(topNo,:)=[1,subid1,subid2];
                elseif iSid==2  % face perp. to y axis  subid1 -> subid_z; subid2 -> subid_x
                    No2fun.subid(topNo,:)=[subid2,1,subid1];
                else            % face perp. to z axis  subid1 -> subid_x; subid2 -> subid_y
                    No2fun.subid(topNo,:)=[subid1,subid2,1];
                end
                topNo=topNo+1;
            end
        end
    end
end
if ~( No2fun.name(topNo)=='' &&  No2fun.name(topNo-1)=='f' )
    error ('Basis number check error!');
end

% domain mode
No2fun.name(topNo:topNo+NBdomain)='d';
for Did=1:mesh.Ndomains
    for subx=2:K+1
        for suby=2:K+1
            for subz=2:K+1
                fun2No.domain(subx,suby,subz,Did)=topNo;
                No2fun.objid(topNo)=Did;
                No2fun.subid(topNo,:)=[subx,suby,subz];   
                topNo=topNo+1;
            end
        end
    end
end
if ~( topNo==Nbasis+1 )
    error ('Basis number check error!');
end

% get domain info
domainInfo=mesh.domains;
Ndomains=mesh.Ndomains;

% ==================== The following do not read mesh any more ===================================
Mvec=
Svec=
MPintegralList=
if strcmp(boundaryType,'second')
    QintegralList=
end

% enummerate basis
for Did=1:Ndomains
    for ix=0:K+1
        for kx=0:K+1
            % cut branch
            if (1)
                continue;
            end
            for iy=0:K+1
                for ky=0:K+1
                    for iz=0:K+1
                        for kz=0:K+1

                
                
                        end
                    end
                end
            end
        end
    end
end


numel(MPintegralList);
end

