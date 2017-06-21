function [ gp_x, gw ] = getGaussPts( obj, ngp )
%Find the gauss points and weigths for Gauss¨CLegendre quadrature
% [George-Gate @2017-06-18]
%
% --------------- Generate gauss points and weights --------------
% gp_x: gauss points; gw: weights

    hisPos=mod(ngp,obj.gHisLen)+1;
    % look in history first
    if length(obj.gw{hisPos})==ngp
        gp_x=obj.gp_x{hisPos};
        gw=obj.gw{hisPos};
    else  % if not found, generate it

        m=(1:ngp-1)';
        c=m./sqrt(4*m.*m-1);
        Jm=spdiags([[c;0],[0;c]],[-1,1],ngp,ngp);
        [Vec,Dig]=eigs(Jm,ceil(ngp/2),'sa');
        gp_x=diag(Dig);
        % normalize V
        Vec=Vec./repmat(sqrt(sum(Vec.*Vec,1)),ngp,1);
        gw=2*abs(Vec(1,:).*Vec(1,:))';
        if mod(ngp,2)==1
            gp_x=[gp_x;-gp_x(end-1:-1:1)];
            gw=  [gw  ; gw(  end-1:-1:1)];
        else
            gp_x=[gp_x;-gp_x(end:-1:1)];
            gw=  [gw  ; gw(  end:-1:1)];
        end
        
        % record to history
        obj.gp_x{hisPos}=gp_x;
        obj.gw{hisPos}=gw;
        
    end
end

