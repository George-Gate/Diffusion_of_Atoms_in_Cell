function plotMesh(obj,displayText)
% Visualize mesh
% Boundary elements and inner elements are marked by different color.
    if ~exist('displayText','var')
        displayText=0;
    end
    
    figure();
    %% plot node
    idOnB=find(obj.nodes.onB == 1);
    idOffB=find(obj.nodes.onB == 0);
    x=obj.nodes.x;
    y=obj.nodes.y;
    z=obj.nodes.z;
    scatter3(x(idOnB ),y(idOnB ),z(idOnB),'c^');hold on;
    scatter3(x(idOffB),y(idOffB),z(idOffB),'bo');

    if displayText
        % Disp the id of node
        text(x,y,z,num2str((1:obj.Nnodes)'),'VerticalAlignment','bottom','color','r');
    end

    %% plot edge
    idOnB=find(obj.edges.onB == 1);
    idOffB=find(obj.edges.onB == 0);
    clear x y z;
    x(1,:)=obj.nodes.x(obj.edges.n(1,:))';y(1,:)=obj.nodes.y(obj.edges.n(1,:))';z(1,:)=obj.nodes.z(obj.edges.n(1,:))';
    x(2,:)=obj.nodes.x(obj.edges.n(2,:))';y(2,:)=obj.nodes.y(obj.edges.n(2,:))';z(2,:)=obj.nodes.z(obj.edges.n(2,:))';
    line(x(:,idOnB ),y(:,idOnB ),z(:,idOnB ),'color','c');
    line(x(:,idOffB),y(:,idOffB),z(:,idOffB),'color','b');

    if displayText
        % Disp the id of edge
        text(sum(x,1)'/2,sum(y,1)'/2,sum(z,1)'/2,num2str((1:obj.Nedges)'),'color','b','HorizontalAlignment','center');
    end



    %% plot surface
    idOnB=find(obj.surfaces.onB == 1);
    idOffB=find(obj.surfaces.onB == 0);
    clear x y z;
    x(1,:)=obj.nodes.x(obj.surfaces.n(1,:))';y(1,:)=obj.nodes.y(obj.surfaces.n(1,:))';z(1,:)=obj.nodes.z(obj.surfaces.n(1,:))';
    x(2,:)=obj.nodes.x(obj.surfaces.n(2,:))';y(2,:)=obj.nodes.y(obj.surfaces.n(2,:))';z(2,:)=obj.nodes.z(obj.surfaces.n(2,:))';
    x(3,:)=obj.nodes.x(obj.surfaces.n(3,:))';y(3,:)=obj.nodes.y(obj.surfaces.n(3,:))';z(3,:)=obj.nodes.z(obj.surfaces.n(3,:))';
    x(4,:)=obj.nodes.x(obj.surfaces.n(4,:))';y(4,:)=obj.nodes.y(obj.surfaces.n(4,:))';z(4,:)=obj.nodes.z(obj.surfaces.n(4,:))';

    fill3(x(:,idOnB),y(:,idOnB),z(:,idOnB),'c','FaceAlpha',0.2,'LineStyle','none')
    fill3(x(:,idOffB),y(:,idOffB),z(:,idOffB),'r','FaceAlpha',0.2,'LineStyle','none')

    if displayText
        x=x([1,3],:);
        y=y([1,3],:);
        z=z([1,3],:);
        text(sum(x,1)'/2,sum(y,1)'/2,sum(z,1)'/2,num2str((1:obj.Nsurfaces)'),'color','k','HorizontalAlignment','center');
    end

    %% plot domain
    if displayText
        x=reshape(obj.domains.xyz(1,:,:),2,obj.Ndomains);
        y=reshape(obj.domains.xyz(2,:,:),2,obj.Ndomains);
        z=reshape(obj.domains.xyz(3,:,:),2,obj.Ndomains);
        text(sum(x,1)'/2,sum(y,1)'/2,sum(z,1)'/2,num2str((1:obj.Ndomains)'),'color','y','HorizontalAlignment','center');
    end


    clear x y z;
    %%
    set(gca,'xTick',[0 1],'xTickMode','manual');
    set(gca,'yTick',[0 1],'yTickMode','manual');
    set(gca,'zTick',[0 1],'zTickMode','manual');
    xlabel('x');ylabel('y');zlabel('z');
end