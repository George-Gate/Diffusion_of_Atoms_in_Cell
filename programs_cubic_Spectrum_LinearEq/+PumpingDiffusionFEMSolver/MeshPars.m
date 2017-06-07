classdef MeshPars
    %MeshPars   Provide parameters for @Mesh3D
    % Example: xMajor=[1,2,5], Nx=[2,3]. Then the x coordinates the mesh points are 1,1.5,2,3,4,5
    properties
        xMajor=[-1,1];   % Should be a column vector with ascending order.
        yMajor=[-1,1];
        zMajor=[-1,1];
        Nx=1;            % Should have length(xMajor)==length(Nx)+1
        Ny=1;
        Nz=1;
    end
    
    properties(Dependent)
        xList
        yList
        zList
        numMeshDomains
    end
    
    methods
        function obj=set.xMajor(obj,val)
            if (isvector(val) && length(val)>=2)
                obj.xMajor=sort(reshape(val,1,length(val)));
                checkDump(obj.xMajor);
            else
                error('xMajor should be a vector with length > 1.');
            end
        end
        function obj=set.yMajor(obj,val)
            if (isvector(val) && length(val)>=2)
                obj.yMajor=sort(reshape(val,1,length(val)));
                checkDump(obj.yMajor);
            else
                error('yMajor should be a vector with length > 1.');
            end
        end
        function obj=set.zMajor(obj,val)
            if (isvector(val) && length(val)>=2)
                obj.zMajor=sort(reshape(val,1,length(val)));
                checkDump(obj.zMajor);
            else
                error('zMajor should be a vector with length > 1.');
            end
        end
        
        function xList=get.xList(obj)
            xList=makeAxisList(obj.xMajor,obj.Nx,'x');
        end
        function yList=get.yList(obj)
            yList=makeAxisList(obj.yMajor,obj.Ny,'y');
        end
        function zList=get.zList(obj)
            zList=makeAxisList(obj.zMajor,obj.Nz,'z');
        end
        function NN=get.numMeshDomains(obj)
            NN=sum(obj.Nx)*sum(obj.Ny)*sum(obj.Nz);
        end
    end
    
end

function checkDump(val)
    % check if there is two number in val that is exactly the same.
    for i=1:length(val)-1
        if (val(i)==val(i+1))
            error('One point should not appear twice!');
        end
    end
end


function coList=makeAxisList(Ps,N,axisName)
    % check the length of majorPoint and N
    if (length(Ps) ~= length(N)+1)
        error(['Input should satisties length(',axisName,'Major) == length(N',axisName,') + 1.']);
    end
    coList=0;
    for i=1:length(N)
        coList=[coList(1:end-1);linspace(Ps(i),Ps(i+1),N(i)+1)'];
    end
end