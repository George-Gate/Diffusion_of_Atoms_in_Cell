function [ vecF ] = projection( obj,fun,maxOrder,ngp )
%projection   Calculate the inner product of base function and a given function
%   vecF is a [maxOrder x 1] column vector with vecF_i=(v_i,fun)
%   fun should be function handle with only one input argument
%   Input argument 'ngp' is optional. It tell the method how many gauss point should be used to calc numerical intergrals.
%   If ngp is not given, the default value is ngp=maxOrder+100.

    if nargin<4
        ngp=maxOrder+100;
    end
    % check the form of fun
    if ~(isa(fun,'function_handle') &&  nargin(fun)==1)
        error('The ''fun'' should be a function handle with exactly one input argument.');
    end
    % get gauss points
    [ gp_x, gw ] = obj.getGaussPts( ngp );

    % evaluate base function at gauss points
    baseVal=obj.funVal(gp_x,1:maxOrder);
    % evaluate fun at gauss points
    try
        fun_num=fun(gp_x);
    catch ME
        fun_num=zeros(size(gp_x));
        for i=1:length(gp_x)
            fun_num(i)=fun(gp_x(i));
        end
    end
    % calc integrals
    fun_num=gw.*fun_num;
    vecF=baseVal'*fun_num;
end

