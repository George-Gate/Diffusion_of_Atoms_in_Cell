figure();
displayText=0;
displayRelation=0;   % this function not finished
%% plot node
idOnB=find(mesh0.nodes.onB == 1);
idOffB=find(mesh0.nodes.onB == 0);
x=mesh0.nodes.x;
y=mesh0.nodes.y;
z=mesh0.nodes.z;
scatter3(x(idOnB ),y(idOnB ),z(idOnB),'c^');hold on;
scatter3(x(idOffB),y(idOffB),z(idOffB),'bo');

if displayText
    % Disp the id of node
    text(x,y,z,num2str((1:mesh0.Nnodes)'),'VerticalAlignment','bottom','color','r');
%     if displayRelation
%         % Disp nodes.s
%         text(x,y,num2str(mesh0.nodes.s(1,:)','\\qquad %i'),'HorizontalAlignment','left','VerticalAlignment','bottom','color','r','interpreter','latex');
%         text(x,y,num2str(mesh0.nodes.s(2,:)','%i \\quad'),'HorizontalAlignment','right','VerticalAlignment','bottom','color','r','interpreter','latex');
%         text(x,y,num2str(mesh0.nodes.s(3,:)','%i \\quad'),'HorizontalAlignment','right','VerticalAlignment','top','color','r','interpreter','latex');
%         text(x,y,num2str(mesh0.nodes.s(4,:)','\\quad %i'),'HorizontalAlignment','left','VerticalAlignment','top','color','r','interpreter','latex');
%     end
end

%% plot edge
idOnB=find(mesh0.edges.onB == 1);
idOffB=find(mesh0.edges.onB == 0);
clear x y z;
x(1,:)=mesh0.nodes.x(mesh0.edges.n(1,:))';y(1,:)=mesh0.nodes.y(mesh0.edges.n(1,:))';z(1,:)=mesh0.nodes.z(mesh0.edges.n(1,:))';
x(2,:)=mesh0.nodes.x(mesh0.edges.n(2,:))';y(2,:)=mesh0.nodes.y(mesh0.edges.n(2,:))';z(2,:)=mesh0.nodes.z(mesh0.edges.n(2,:))';
line(x(:,idOnB ),y(:,idOnB ),z(:,idOnB ),'color','c');
line(x(:,idOffB),y(:,idOffB),z(:,idOffB),'color','b');

if displayText
    if displayRelation
        % Disp the id of edge and adjacent surfaces
        text(sum(x,1)'/2,sum(y,1)'/2,num2str([mesh0.edges.s(1,:)',(1:mesh0.Nedges)',mesh0.edges.s(2,:)']),'color','b','HorizontalAlignment','center');
    else
        % Disp the id of edge
        text(sum(x,1)'/2,sum(y,1)'/2,sum(z,1)'/2,num2str((1:mesh0.Nedges)'),'color','b','HorizontalAlignment','center');
    end
end



%% plot surface
idOnB=find(mesh0.surfaces.onB == 1);
idOffB=find(mesh0.surfaces.onB == 0);
clear x y z;
x(1,:)=mesh0.nodes.x(mesh0.surfaces.n(1,:))';y(1,:)=mesh0.nodes.y(mesh0.surfaces.n(1,:))';z(1,:)=mesh0.nodes.z(mesh0.surfaces.n(1,:))';
x(2,:)=mesh0.nodes.x(mesh0.surfaces.n(2,:))';y(2,:)=mesh0.nodes.y(mesh0.surfaces.n(2,:))';z(2,:)=mesh0.nodes.z(mesh0.surfaces.n(2,:))';
x(3,:)=mesh0.nodes.x(mesh0.surfaces.n(3,:))';y(3,:)=mesh0.nodes.y(mesh0.surfaces.n(3,:))';z(3,:)=mesh0.nodes.z(mesh0.surfaces.n(3,:))';
x(4,:)=mesh0.nodes.x(mesh0.surfaces.n(4,:))';y(4,:)=mesh0.nodes.y(mesh0.surfaces.n(4,:))';z(4,:)=mesh0.nodes.z(mesh0.surfaces.n(4,:))';

fill3(x(:,idOnB),y(:,idOnB),z(:,idOnB),'c','FaceAlpha',0.2,'LineStyle','none')
fill3(x(:,idOffB),y(:,idOffB),z(:,idOffB),'r','FaceAlpha',0.2,'LineStyle','none')

if displayText
    if displayRelation
        for i=1:mesh0.Nsurfaces
            x=sum(mesh0.surfaces.x(:,i))/2;
            y=sum(mesh0.surfaces.y(:,i))/2;
            text(x,y,num2str([mesh0.surfaces.n(2,i),mesh0.surfaces.e(1,i),mesh0.surfaces.n(1,i);  ...
                              mesh0.surfaces.e(2,i),        i            ,mesh0.surfaces.e(4,i);       ...
                              mesh0.surfaces.n(3,i),mesh0.surfaces.e(3,i),mesh0.surfaces.n(4,i)]),'color','k','HorizontalAlignment','center');
        end
    else
        x=x([1,3],:);
        y=y([1,3],:);
        z=z([1,3],:);
        text(sum(x,1)'/2,sum(y,1)'/2,sum(z,1)'/2,num2str((1:mesh0.Nsurfaces)'),'color','k','HorizontalAlignment','center');
    end
end

%% plot domain
if displayText
    x=reshape(mesh0.domains.xyz(1,:,:),2,mesh0.Ndomains);
    y=reshape(mesh0.domains.xyz(2,:,:),2,mesh0.Ndomains);
    z=reshape(mesh0.domains.xyz(3,:,:),2,mesh0.Ndomains);
    text(sum(x,1)'/2,sum(y,1)'/2,sum(z,1)'/2,num2str((1:mesh0.Ndomains)'),'color','y','HorizontalAlignment','center');
end


clear x y z;
%%
set(gca,'xTick',[0 1],'xTickMode','manual');
set(gca,'yTick',[0 1],'yTickMode','manual');
set(gca,'zTick',[0 1],'zTickMode','manual');
xlabel('x');ylabel('y');zlabel('z');