% fun_fem1.m 计算有限元解, 使用线性基函数空间
function u_fem = fun_fem1(N,fem_obj)
% 使用 fem_class 类获取边值问题参数
a=fem_obj.a; %左端点
b=fem_obj.b; %右端点

%% 参数初始化设置
h=(b-a)/N; %步长
nodes=a:h:b; %节点

%初始化刚度矩阵和载荷向量
n=N; %基函数数量
K=zeros(n,n); %刚度矩阵
F=zeros(n,1); %载荷向量

%% 单元刚度矩阵和单元载荷向量的生成与组装
for i=1:N
    % 当前单元e_i: [x_{i-1}, x_i]
    node_left = nodes(i);
    % node_right = nodes(i+1);

    % 确定当前单元相关的基函数索引
    if i==1
        % 第一个单元 [x_0, x_1]: 只涉及基函数 \phi_1
        indices = 1;
        num_local_bases = 1;
    else
        % 中间单元 [x_{i-1}, x_i]: 涉及基函数 \phi_{i-1} 和 \phi_i
        indices = [i-1, i];
        num_local_bases = 2;
    end

    % 单元刚度矩阵(随基函数选取变化)
    if num_local_bases == 1
        Ke=1/h;
    else
        Ke=[1/h -1/h; -1/h 1/h];
    end

    % 单元载荷向量
    %Gauss 积分点和权重
    gauss_points = [node_left + (1-1/sqrt(3))*h/2; node_left + (1+1/sqrt(3))*h/2];
    gauss_weights = [h/2, h/2];

    %计算基函数在高斯点的值(随基函数选取变化)
    if num_local_bases ==1
        phi_val=[(1-1/sqrt(3))/2; (1+1/sqrt(3))/2]; % φ_1 在两个高斯点的值
    else
        phi_val_1=[(1+1/sqrt(3))/2; (1-1/sqrt(3))/2]; % φ_{i-1} 在两个高斯点的值
        phi_val_2=[(1-1/sqrt(3))/2; (1+1/sqrt(3))/2]; % φ_i 在两个高斯点的值
    end

    % 组装单元载荷向量
    if num_local_bases == 1
        % 第一个单元 [x_0, x_1]: 只涉及基函数 \phi_1
        Fe=gauss_weights*(fem_obj.fun_source(gauss_points).*phi_val);
    else
        % 中间单元 [x_{i-1}, x_i]: 涉及基函数 \phi_{i-1} 和 \phi_i
        Fe_1=gauss_weights*(fem_obj.fun_source(gauss_points).*phi_val_1);
        Fe_2=gauss_weights*(fem_obj.fun_source(gauss_points).*phi_val_2);
        Fe=[Fe_1; Fe_2];
    end

    % 将单元刚度矩阵和载荷向量组装到全局矩阵和向量中
    K(indices, indices) = K(indices, indices) + Ke;
    F(indices) = F(indices) + Fe;
end

%求解线性方程组
u_coeff = K \ F;

% 构建完整的解向量 (n+1 个点)
u_fem = [0; u_coeff];
end