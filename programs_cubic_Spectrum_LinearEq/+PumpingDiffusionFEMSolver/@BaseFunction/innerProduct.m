function [ M ] = innerProduct( obj,weightFun,maxOrder,ngp )
%innerProduct  Calculate weighted inner product between basis functions
%   M is a [maxOrder x maxOrder] matrix with M_ik=(v_i,weightFun*v_k)
%   weightFun should be function handles with only one input argument
%   Input argument 'ngp' is optional. It tell the method how many gauss point should be used to calc numerical intergrals.
%   If ngp is not given, the default value is ngp=2*maxOrder+100 .
    if nargin<4
        ngp=2*maxOrder+100;
    end
    % check the form of weightFun
    if ~(isa(weightFun,'function_handle') &&  nargin(weightFun)==1)
        error('The ''weightFun'' should be a function handle with exactly one input argument.');
    end
    % get gauss points
    [ gp_x, gw ] = obj.getGaussPts( ngp );
    % evaluate base function at gauss points
    baseVal=obj.funVal(gp_x,1:maxOrder);
    % evaluate weightFun at gauss points
    try
        wFun_num=weightFun(gp_x);
    catch ME
        wFun_num=zeros(size(gp_x));
        for i=1:length(gp_x)
            wFun_num(i)=weightFun(gp_x(i));
        end
    end
    % calc integrals
    wFun_num=gw.*wFun_num;
    wFun_num=repmat(wFun_num,1,maxOrder);
    M=baseVal'*(wFun_num.*baseVal);
end

