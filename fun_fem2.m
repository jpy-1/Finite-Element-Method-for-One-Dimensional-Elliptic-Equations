% fun_fem2.m 计算有限元解，使用二次基函数空间
function u_fem=fun_fem2(N,fem_obj)
%% - 边值问题：u''=f 
% 使用 fem_class 类获取边值问题参数   
a=fem_obj.a; %左端点
b=fem_obj.b; %右端点
%% 参数初始化设置
h=(b-a)/N; %步长
nodes=zeros(2*N+1,1); %节点
for i=1:N
    nodes(2*i-1)=a+(i-1)*h; %偶数节点
    nodes(2*i)=a+(i-1)*h+h/2; %奇数节点
end
nodes(2*N+1)=b; %最后一个节点
n=2*N; % 基函数数量

%初始化刚度矩阵和载荷向量
K=zeros(n,n); %刚度矩阵
F=zeros(n,1); %载荷向量

%% 单元刚度矩阵和单元载荷向量的生成与组装
for i=1:N
    node_left = nodes(2*i-1);
    % node_mid = nodes(2*i);
    % node_right = nodes(2*i+1);
    % 确定当前单元相关的基函数索引
    if i==1
        % 第一个单元 [x_0, x_1]: 只涉及基函数 \phi_1
        indices = [1,2];
        num_local_bases = 2;
    else
        % 中间单元 [x_{i-1}, x_i]: 涉及基函数 \phi_{i-1} 和 \phi_i
        indices = [2*i-2, 2*i-1, 2*i];
        num_local_bases = 3;
    end

    if num_local_bases == 2
        % 单元刚度矩阵
        Ke=[16/(3*h),-8/(3*h);
        -8/(3*h),7/(3*h)];
    else
        % 单元刚度矩阵
        Ke=[7/(3*h) -8/(3*h) 1/(3*h);
            -8/(3*h) 16/(3*h) -8/(3*h);
            1/(3*h) -8/(3*h) 7/(3*h)];
    end

    % 单元载荷向量
    %Gauss 积分点和权重
    gauss_points = [node_left + (1-1/sqrt(3))*h/2; node_left + (1+1/sqrt(3))*h/2];
    gauss_weights = [h/2, h/2];
    if num_local_bases ==2
        phi_val_1=[2/3; 2/3]; % φ_{i-1/2} 在两个高斯点的值
        phi_val_2=[(1-sqrt(3))/6; (1+sqrt(3))/6]; % φ_i 在两个高斯点的值
    else
        phi_val_1=[(sqrt(3)+1)/6; (1-sqrt(3))/6]; % φ_{i-1} 在两个高斯点的值
        phi_val_2=[2/3; 2/3]; % φ_i-1/2 在两个高斯点的值
        phi_val_3=[(1-sqrt(3))/6; (1+sqrt(3))/6]; % φ_{i} 在两个高斯点的值
    end

    % 组装单元载荷向量
    if num_local_bases == 2
        % 第一个单元 [x_0, x_1]: 只涉及基函数 \phi_1
        Fe_1=gauss_weights*(fem_obj.fun_source(gauss_points).*phi_val_1);
        Fe_2=gauss_weights*(fem_obj.fun_source(gauss_points).*phi_val_2);
        Fe=[Fe_1; Fe_2];
    else
        % 中间单元 [x_{i-1}, x_i]: 涉及基函数 \phi_{i-1} 和 \phi_i
        Fe_1=gauss_weights*(fem_obj.fun_source(gauss_points).*phi_val_1);
        Fe_2=gauss_weights*(fem_obj.fun_source(gauss_points).*phi_val_2);
        Fe_3=gauss_weights*(fem_obj.fun_source(gauss_points).*phi_val_3);
        Fe=[Fe_1; Fe_2; Fe_3];
    end

    % 将单元刚度矩阵和载荷向量组装到全局矩阵和向量中
    K(indices, indices) = K(indices, indices) + Ke;
    F(indices) = F(indices) + Fe;
end

%求解线性方程组
u_coeff = K \ F;

u_fem=[0;u_coeff]; %构建完整的解向量 (2n+1 个点)

end
