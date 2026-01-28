% fem_class.m 封装边值问题的基本属性，其中u(a)=0, u'(b)=0：
% 区间端点 a, b
% 源项函数句柄 f(x)
% 解析解函数句柄 u(x)

 classdef fem_class

    properties (Constant)  % 声明为常量属性（不变的量）
        a = 0;             % 左端点（固定值）
        b = pi/8;          % 右端点（固定值）
    end
    
    methods
        % 解析解函数：u(x) = sin(4x)/16（问题解析解）
        function y = fun_analytical(~, x)
            y = sin(4 .* x) / 16;
        end
        
        % 源函数：f(x) = sin(4x)（问题源项）
        function f = fun_source(~, x)
            f = sin(4 .* x);
        end
    end
end