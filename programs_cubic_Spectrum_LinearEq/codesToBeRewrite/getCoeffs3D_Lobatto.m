function [ Mvec, Svec, fun2IntegralID_MP, MPid, Nbasis, fun2No, No2fun, getNoByIxyz ] = getCoeffs3D_Lobatto( boundaryType, mesh, K )
% Generate the coeff matrices for a box area.
%    Do not support non constant first kind boundary condition so far. 
% [Inputs]
%    boundaryType: 'first', 'second', 'firstUniform', 'secondZero'
%    if boundaryType == 'firstUniform', we will give a matrix that solves zero boundary problem. i.e. our solution is (u-c) where c is the boundary value.
%
%    K: the cutoff of the order of Lobatto polynomials
%
% [Outputs]
%   fun2IntegralID_MP: see below
%   MPid: a list of integral to be calculated. matrix MP depends on these integrals.
%   Mvec, Svec
%
% 
% #TODO: faceDir, edgeDir currently not used, delete them?


% ==================== read mesh ================================
if strcmp(boundaryType,'firstUniform') || strcmp(boundaryType,'first')
    boundaryBasis=0;
elseif strcmp(boundaryType,'secondZero') || strcmp(boundaryType,'second')
    boundaryBasis=1;
else
    error(['Unknow boundary type: ',boundaryType]);
    boundaryBasis=0;
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
fun2No.node(No2fun.objid(1:NBnode))=(1:NBnode)';

% edge mode
No2fun.name(NBnode+1:NBnode+NBedge)='e';
topNo=NBnode+1;
for Nid=1:mesh.Nnodes
    for iEid=1:3
        Eid=mesh.nodes.e(iEid,Nid);
        if (Eid==0); continue; end
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
if topNo>1 && ~( No2fun.name(topNo)==0 &&  No2fun.name(topNo-1)=='e' )
    error ('Basis number check error!');
end

% face mode
No2fun.name(topNo:topNo+NBface-1)='f';
for Nid=1:mesh.Nnodes
    for iSid=1:3
        Sid=mesh.nodes.s(iSid,Nid);
        if (Sid==0); continue; end
        faceDir(Sid)=iSid;
        if (mesh.surfaces.onB(Sid) && ~boundaryBasis)
            % skip boundary elements if don't need boundary basis
            continue;
        end
        for subid1=2:K+1
            for subid2=2:K+1
                fun2No.face(subid1,subid2,Sid)=topNo;
                No2fun.objid(topNo)=Sid;
                if iSid==1      % face perp. to x axis  subid1 -> subid_y; subid2 -> subid_z
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
if topNo>1 && ~( No2fun.name(topNo)==0 &&  No2fun.name(topNo-1)=='f' )
    error ('Basis number check error!');
end

% domain mode
No2fun.name(topNo:topNo+NBdomain-1)='d';
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
% get boundary elements
if boundaryBasis
    NonB=mesh.nodes.onB;
    EonB=mesh.edges.onB;
    SonB=mesh.surfaces.onB;
end

% ==================== The following do not read mesh any more ===================================
% ---------------------------- Create getNoByIxyz for heavily speed up ---------------------------------
getNoByIxyz=zeros(K+2,K+2,K+2,Ndomains);
for Did=1:Ndomains
    for iz=0:K+1                % x basis for row
        for iy=0:K+1            % y basis for row
            for ix=0:K+1        % z basis for row
                getNoByIxyz(ix+1,iy+1,iz+1,Did)=getNoByIxyz_fun(ix,iy,iz,Did,domainInfo,fun2No);
            end
        end
    end
end
% ---------------------------- Generate Mvec ans Svec --------------------------------------------
Mvec=zeros(Ndomains*(3*K+8)^3+1,3);
Svec=zeros(Ndomains*(3*K+8)^3+1,3);
topM=1; topS=1;

% enummerate basis to get Mvec and Svec
for Did=1:Ndomains
    h=domainInfo.h(1:3,Did);
    for ix=0:K+1           % x basis for row
        for kx=max(ix-3,0):min(ix+3,K+1)   % x basis for column
            % cut branch
            if zeroResult(ix,kx)
                continue;
            end
%             disp([num2str([ix,kx])]);
            phiphi_x=phiphi(ix,kx,h(1));
            dphidphi_x=dphidphi(ix,kx,h(1));
            for iy=0:K+1  % y basis for row
                for ky=max(iy-3,0):min(iy+3,K+1)
                    % cut branch
                    if zeroResult(iy,ky)
                        continue;
                    end
                    phiphi_y=phiphi(iy,ky,h(2));
                    dphidphi_y=dphidphi(iy,ky,h(2));
                    for iz=0:K+1    % z basis for row
%                         rowNo=getNoByIxyz_fun(ix,iy,iz,Did,domainInfo,fun2No);
                        rowNo=getNoByIxyz(ix+1,iy+1,iz+1,Did);
                        if (rowNo==0)
                            continue;
                        end
                        for kz=max(iz-3,0):min(iz+3,K+1)
                            % cut branch
                            if zeroResult(iz,kz)
                                continue;
                            end
%                             colNo=getNoByIxyz_fun(kx,ky,kz,Did,domainInfo,fun2No);
                            colNo=getNoByIxyz(kx+1,ky+1,kz+1,Did);
                            if (colNo==0)
                                continue;
                            end
                            phiphi_z=phiphi(iz,kz,h(3));
                            dphidphi_z=dphidphi(iz,kz,h(3));
                            % set M matrix
                            Mvec(topM,1)=rowNo;
                            Mvec(topM,2)=colNo;
                            Mvec(topM,3)=phiphi_x*phiphi_y*phiphi_z;
                            if abs(Mvec(topM,3))>0
                                topM=topM+1;
                            end
                            % set S matrix
                            Svec(topS,1)=rowNo;
                            Svec(topS,2)=colNo;
                            Svec(topS,3)=dphidphi_x*phiphi_y*phiphi_z ...
                                        +phiphi_x*dphidphi_y*phiphi_z ...
                                        +phiphi_x*phiphi_y*dphidphi_z;
                            if abs(Svec(topS,3))>0
                                topS=topS+1;
                            end
                        end
                    end
                end
            end
        end
    end
end

% trim Mvec & Svec
Mvec=Mvec(1:topM-1,1:3);
Svec=Svec(1:topS-1,1:3);

% -------------------------- Generate integral list for calculating MP and Q --------------------------
% assume that the weight function P(r) is separable as P(r)=A(x)B(y)C(z)
% generate single axis integral list for MP
fun2IntegralID_MP=reshape(1:(K+2)*(K+2)*3*Ndomains,K+2,K+2,3,Ndomains);
% fun2IntegralID(l,m,k,Did) is the No. of integral (phi_l,A(x)*phi_m) on domain Did if k=1. k=2/3 refer to y/z axis integral.


MPid=zeros(Ndomains*(K+2)^6,5);   % MPid(:,1/2): rowNo/colNo; MPid(:,3/4/5): x/y/z part of integral. 
topMP=1;

% enummerate basis to get MPid
for Did=1:Ndomains
    for iz=0:K+1                % x basis for row
        for iy=0:K+1            % y basis for row
            for ix=0:K+1        % z basis for row
%                 rowNo=getNoByIxyz_fun(ix,iy,iz,Did,domainInfo,fun2No);
                rowNo=getNoByIxyz(ix+1,iy+1,iz+1,Did);
                if rowNo==0 ; continue; end
                for kz=0:K+1    % x basis for column
                    for ky=0:K+1
                        for kx=0:K+1
%                             colNo=getNoByIxyz_fun(kx,ky,kz,Did,domainInfo,fun2No);
                            colNo=getNoByIxyz(kx+1,ky+1,kz+1,Did);
                            if colNo==0 ; continue; end
                            MPid(topMP,1)=rowNo;
                            MPid(topMP,2)=colNo;
                            MPid(topMP,3)=fun2IntegralID_MP(ix+1,kx+1,1,Did);
                            MPid(topMP,4)=fun2IntegralID_MP(iy+1,ky+1,2,Did);
                            MPid(topMP,5)=fun2IntegralID_MP(iz+1,kz+1,3,Did);
                            topMP=topMP+1;
                        end
                    end
                end
            end
        end
    end
end

% trim MPid
MPid=MPid(1:topMP-1,:);

end

% ============================ Private Functions =======================================

% given subid_x/y/z and domain id, return the relevant basis No.
function No=getNoByIxyz_fun(subx,suby,subz,Did,domainInfo,fun2No)
    No=0;
    type=sum([subx<=1;suby<=1;subz<=1]);
    if (type==3)  % node basis
        Nid=domainInfo.n(subx+1,suby+1,subz+1,Did);
        No=fun2No.node(Nid);
    elseif (type==2)  % edge basis
        if (subx>1)
            Eid=domainInfo.e(suby+1,subz+1,1,Did);
            No=fun2No.edge(subx,Eid);
        elseif (suby>1)
            Eid=domainInfo.e(subz+1,subx+1,2,Did);
            No=fun2No.edge(suby,Eid);
        else  % subz>1
            Eid=domainInfo.e(subx+1,suby+1,3,Did);
            No=fun2No.edge(subz,Eid);
        end
    elseif (type==1)  % face basis
        if (subx<=1)
            Sid=domainInfo.s(subx+1,1,Did);
            No=fun2No.face(suby,subz,Sid);
        elseif (suby<=1)
            Sid=domainInfo.s(suby+1,2,Did);
            No=fun2No.face(subz,subx,Sid);
        else   % subz<=1
            Sid=domainInfo.s(subz+1,3,Did);
            No=fun2No.face(subx,suby,Sid);
        end
    else             % domain basis
        No=fun2No.domain(subx,suby,subz,Did);
    end
end

% return integral (phi_l,phi_m), where li/mk=0 stands for phi_1(x;a); li/mk=1 stands for phi_m(x;b)
function v=phiphi(li,mk,h)
    v=0;
    if (li<0 || mk<0); return; end
    if (li>1 && mk>1)
        if mk==li
            v=h/(2*li-1)/(2*li-1)*(1/(2*li+1)+1/(2*li-3));
        elseif mk==li+2
            v=-h/(2*li-1)/(2*li+1)/(2*li+3);
        elseif mk==li-2
            v=-h/(2*li-1)/(2*li-3)/(2*li-5);
        end
    elseif (li>1)  % mk<=1
        if (mk==0 && li<=3)
            v=h/(36*li-78);    % -h/6 for l=2 and h/30 for l=3.
        elseif (mk==1 && li<=3)
            v=-h/(24*li-42);    % -h/6 for l=2 and -h/30 for l=3.
        end
    elseif (mk>1)  % li<=1
        if (li==0 && mk<=3)
            v=h/(36*mk-78);    % -h/6 for m=2 and h/30 for m=3.
        elseif (li==1 && mk<=3)
            v=-h/(24*mk-42);    % -h/6 for m=2 and -h/30 for m=3.
        end
    else           % 0<=li<=1 && 0<=mk<=1
        v=h/3/(abs(li-mk)+1);   % h/3 for li==mk; h/6 for li~=mk
    end
end

% return integral (dphi_l,dphi_m), where li/mk=0 stands for phi_1(x;a); li/mk=1 stands for phi_m(x;b)
function v=dphidphi(li,mk,h)
    v=0;
    if (li<0 || mk<0); return; end
    if (li>1 && mk>1)
        if (li==mk)
            v=4/h/(2*li-1);
        end
    elseif (li<=1 && mk<=1)
        v=(1-2*abs(li-mk))/h;
    end

end

% decide whether a conbination of (l,m) give non-zero integral (phi_l,phi_m)
% used to cut branch in basis enumeration
function con=zeroResult(l,m)
    % l>=4
    if (l>=4 && ~(m==l-2 || m==l || m==l+2))
        con=1;return;
    end
    % l==1
    if (l==1 && m==4)
        con=1;return;
    end
    % l=2/3
    if (l==2 || l==3) && ~(m==0 || m==1 || m==l || m==l+2)
        con=1;return;
    end
    con=0;
end