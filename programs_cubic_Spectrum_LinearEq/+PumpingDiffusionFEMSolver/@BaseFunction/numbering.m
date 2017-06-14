function Nbasis=numbering( obj,mesh,maxOrder )
%numbering  Generate the global No. of each base for a given @Mesh3D object.
%    maxOrder should be at least 2.
%    Assume that v_1(-1)~=0, v_1(1)==0;  v_2(-1)==0, v_2(1)~=0;  v_i(+-1)==0 for i>2
%    Without this assumption, the numbering is wrong, consider overload this method in that case.

    if maxOrder<2
        error('maxOrder should be at least 2.');
    end
    
    startT=tic;
    
    K=maxOrder-2;
% ----------------------- Set fun2No and No2fun -------------------------------------------
    NBnode=mesh.Nnodes;
    NBedge=mesh.Nedges*K;
    NBface=mesh.Nsurfaces*K*K;
    NBdomain=mesh.Ndomains*K*K*K;
    Nbasis=NBnode+NBedge+NBface+NBdomain;

    % array def.
    fun2No.node=zeros(mesh.Nnodes,1);
    fun2No.edge=zeros(K+2,mesh.Nedges);              % subid=1,2: base that are non-zero on x=+-1 ;  subid=3~K+2: base that vanish at boundary.  
                                                     % fun2No.edge(1/2,Eid) is never refered.
    fun2No.face=zeros(K+2,K+2,mesh.Nsurfaces);       % fun2No.face(subid_x,subid_y,Fid)/(subid_y,subid_z,Fid)/(subid_z,subid_x,Fid) for different face direction
    fun2No.domain=zeros(K+2,K+2,K+2,mesh.Ndomains);  % fun2No.domain(subid_x, subid_y, subid_z, Did)
    No2fun.name=char(zeros(Nbasis,1));   % possible values: 'n' for node; 'e' for edge; 'f' for face; 'd' for domain
    No2fun.objid=zeros(Nbasis,1);
    No2fun.subid=zeros(Nbasis,3);        % subid=[subid_x;subid_y;subid_z], subid=0 means that that subid is meaningless for the current basis.

    % node mode
    No2fun.name(1:NBnode)='n';
    No2fun.subid(1:NBnode,1:3)=0;   % subid has no meaning for node mode basis.
    No2fun.objid(1:NBnode)=(1:mesh.Nnodes)';
    fun2No.node(No2fun.objid(1:NBnode))=(1:NBnode)';

    % edge mode
    No2fun.name(NBnode+1:NBnode+NBedge)='e';
    topNo=NBnode+1;
    for Nid=1:mesh.Nnodes
        for iEid=1:3
            Eid=mesh.nodes.e(iEid,Nid);
            if (Eid==0); continue; end
            for subid=3:K+2
                fun2No.edge(subid,Eid)=topNo;
                No2fun.objid(topNo)=Eid;
                No2fun.subid(topNo,:)=0;  % subid has no meaning on the direction that perpendicular to the edge
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
            for subid1=3:K+2
                for subid2=3:K+2
                    fun2No.face(subid1,subid2,Sid)=topNo;
                    No2fun.objid(topNo)=Sid;
                    if iSid==1      % face perp. to x axis  subid1 -> subid_y; subid2 -> subid_z
                        No2fun.subid(topNo,:)=[0,subid1,subid2];
                    elseif iSid==2  % face perp. to y axis  subid1 -> subid_z; subid2 -> subid_x
                        No2fun.subid(topNo,:)=[subid2,0,subid1];
                    else            % face perp. to z axis  subid1 -> subid_x; subid2 -> subid_y
                        No2fun.subid(topNo,:)=[subid1,subid2,0];
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
        for subx=3:K+2
            for suby=3:K+2
                for subz=3:K+2
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

% ---------------------------- Create getNoByIxyz ---------------------------------
    getNoByIxyz=zeros(K+2,K+2,K+2,mesh.Ndomains);
    for Did=1:mesh.Ndomains
        for iz=1:K+2                
            for iy=1:K+2           
                for ix=1:K+2      
                    getNoByIxyz(ix,iy,iz,Did)=getNoByIxyz_fun(ix,iy,iz,Did,mesh.domains,fun2No);
                end
            end
        end
    end
% ---------------------------- Store numbering result --------------------------------
    obj.getNoByIxyz=getNoByIxyz;
    obj.fun2No=fun2No;
    obj.No2fun=No2fun;
    obj.meshPars=mesh.meshPars_current;
    obj.maxOrder=maxOrder;    % set maxOrder property to allow other methods to create cache for speed up.
    obj.Nbasis=Nbasis;
    
    import PumpingDiffusionFEMSolver.library.sec2hms;
    disp(['Time used to generate base numbering: ',sec2hms(toc(startT))]);
    
end


% ============================ Private Functions =======================================

% given subid_x/y/z and domain id, return the relevant basis No.
function No=getNoByIxyz_fun(subx,suby,subz,Did,domainInfo,fun2No)
    No=0;
    type=sum([subx<=2;suby<=2;subz<=2]);
    if (type==3)  % node basis
        Nid=domainInfo.n(subx,suby,subz,Did);
        No=fun2No.node(Nid);
    elseif (type==2)  % edge basis
        if (subx>2)
            Eid=domainInfo.e(suby,subz,1,Did);
            No=fun2No.edge(subx,Eid);
        elseif (suby>2)
            Eid=domainInfo.e(subz,subx,2,Did);
            No=fun2No.edge(suby,Eid);
        else  % subz>2
            Eid=domainInfo.e(subx,suby,3,Did);
            No=fun2No.edge(subz,Eid);
        end
    elseif (type==1)  % face basis
        if (subx<=2)
            Sid=domainInfo.s(subx,1,Did);
            No=fun2No.face(suby,subz,Sid);
        elseif (suby<=2)
            Sid=domainInfo.s(suby,2,Did);
            No=fun2No.face(subz,subx,Sid);
        else   % subz<=2
            Sid=domainInfo.s(subz,3,Did);
            No=fun2No.face(subx,suby,Sid);
        end
    else             % domain basis
        No=fun2No.domain(subx,suby,subz,Did);
    end
end

