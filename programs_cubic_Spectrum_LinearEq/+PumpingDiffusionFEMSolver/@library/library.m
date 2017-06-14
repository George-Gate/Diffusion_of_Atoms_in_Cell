classdef library
    %library  A collection of library functions 
    
    properties
    end
    
    methods(Static)
        [ gp_x, gw ] = getGaussPts( ngp );
        [ varargout ] = matStat( A, name );
        I = qgauss( func, x1, x2, ngp );
        hmsStr = sec2hms( tSec );
        mat = vec2full( iList, jList, vList, m, n );
    end
    
end

