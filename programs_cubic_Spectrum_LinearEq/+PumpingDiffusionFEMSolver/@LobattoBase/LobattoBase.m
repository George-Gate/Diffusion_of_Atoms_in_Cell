classdef LobattoBase < PumpingDiffusionFEMSolver.BaseFunction
    %LobattoBase  The base function that use Lobatto function and piecewise linear function.
    %   All basis is defined on [-1,1].
    %   v_1(x)=(1-x)/2;    v_2(x)=(1+x)/2;
    %   v_k(x)=Lo_{k-1}(x),   k>2
    
    
    properties(Access=private)
        phiphiMatrix=[];
        phidphiMatrix=[];
        dphidphiMatrix=[];
    end
    
    methods
        function obj=LobattoBase()
            obj=obj@PumpingDiffusionFEMSolver.BaseFunction();
            obj.basisName='Lobatto';
        end
        
        % return inner product (v_i,v_k)
        function val=phiphi(obj,i,k,h)
            maxik=max(max(i(:)),max(k(:)));
            if ( size(obj.phiphiMatrix,1)<maxik )
                % regenerate phiphiMatrix
                K=max(obj.maxOrder,maxik)+10;  % set the matrix a bit larger
                l=(2:K-1)';
                diagEl=[1/3;1/3;  1./(2*l-1)./(2*l-1).*( 1./(2*l+1) + 1./(2*l-3) )];
                rowID=(1:K)';colID=(1:K)';
                l=(2:K-3)';
                subsubDiag=[-1/6;-1/30;  -1./(2*l-1)./(2*l+1)./(2*l+3)];
                rowID=[rowID;(1:K-2)';(3:K)';1;1;2;2;3;4];
                colID=[colID;(3:K)';(1:K-2)';2;4;1;3;2;1];
                obj.phiphiMatrix=sparse(rowID,colID,[diagEl;subsubDiag;subsubDiag;...
                                 1/6;1/30;1/6;-1/6;-1/6;1/30],K,K);
            end
            val=h*obj.phiphiMatrix(i,k);
        end
        
        % return inner product (v_i',v_k')
        function val=dphidphi(obj,i,k,h)
            maxik=max(max(i(:)),max(k(:)));
            if (size(obj.dphidphiMatrix,1)<maxik)
                % regenerate dphidphiMatrix
                K=max(obj.maxOrder,maxik)+10;  % set the matrix a bit larger
                l=(2:K-1)';
                diagEl=[1;1;4./(2*l-1)];
                rowID=[(1:K)';1;2];
                colID=[(1:K)';2;1];
                obj.dphidphiMatrix=sparse(rowID,colID,[diagEl;-1;-1],K,K);
            end
            val=obj.dphidphiMatrix(i,k)/h;
        end
        
        % return inner product (v_i,v_k')
        function val=phidphi(obj,i,k)
            maxik=max(max(i(:)),max(k(:)));
            if (size(obj.phidphiMatrix,1)<maxik)
                % regenerate dphidphiMatrix
                K=max(obj.maxOrder,maxik)+10;  % set the matrix a bit larger
                l=(2:K-2)';
                subdiag=[1/2;1/3;  2./(2*l-1)./(2*l+1)  ];
                rowID=[(1:K-1)';(2:K)';1;2;1;3];
                colID=[(2:K)';(1:K-1)';1;2;3;1];
                obj.phidphiMatrix=sparse(rowID,colID,[subdiag;-subdiag;-1/2;1/2;-1/3;1/3],K,K);
            end
            val=obj.phidphiMatrix(i,k);
        end
    end
    
    methods(Static)
        % return the value of base function at given points
        function val=funVal(xList,i)
            import PumpingDiffusionFEMSolver.LobattoBase.*
            if ~( isvector(xList) && isvector(i) )
                error('xList and i should be vectors.');
            end
            szX=isrow(xList);
            % make sure xList is a row vector
            if ~szX
                xList=xList';
            end
            j=i-2;
            j(j<1)=1;
            val=lobattoP_N(j,xList)';   % get a result as if xList is a column vector
            for iii=1:length(i)
                if i(iii)==1
                    val(:,iii)=(1-xList')/2;
                elseif i(iii)==2
                    val(:,iii)=(1+xList')/2;
                end
            end
            
            % reshape val
            if szX
                val=val';
            end
        end
        
        % return the value of the first derivative of base function at given points
        function val=funFirstDerivative(xList,i)
            import PumpingDiffusionFEMSolver.LobattoBase.*
            if ~( isvector(xList) && isvector(i) )
                error('xList and i should be vectors.');
            end
            szX=isrow(xList);
            % make sure xList is a row vector
            if ~szX
                xList=xList';
            end
            j=i-2;
            j(j<1)=1;
            val=legendreP_N(j,xList)';   % get a result as if xList is a column vector
            val(:,i==1)=-1/2;
            val(:,i==2)=1/2;
            
            % reshape val
            if szX
                val=val';
            end
        end
        
    end
    
    methods(Static,Hidden,Access=private)
        [ result ] = legendreP_N( nList,xList );
        [ result ] = lobattoP_N( jList,xList );
    end
    
end

