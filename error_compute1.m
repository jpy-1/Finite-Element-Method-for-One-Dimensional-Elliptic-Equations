% error_compute1.m 计算使用线性基函数空间的有限元解的L1, L2, Linf误差
function error = error_compute1(u_fem,fem_obj,N)
    % 初始化误差
    L1 = 0;
    L2 = 0;
    Linf = 0;
    a=fem_obj.a; %左端点
    b=fem_obj.b; %右端点
    h=(b-a)/N; %步长
    nodes=a:h:b; %节点
    sample_num=50;%L^\infty误差单元内采样点数量

    %遍历每个单元计算误差
    for i=1:N
        % 当前单元e_i: [x_{i-1}, x_i]
        node_left = nodes(i);%单元左端点
        node_right = nodes(i+1);%单元右端点
        u_fem_left = u_fem(i);       % 单元左端点的有限元解
        u_fem_right = u_fem(i+1);    % 单元右端点的有限元解

        %使用Gauss积分计算L1,L2误差
        gauss_points = [node_left + (1-1/sqrt(3))*h/2; node_left + (1+1/sqrt(3))*h/2];
        gauss_weights = [h/2, h/2];

        % 计算解析解在高斯点的值
        u_analytical = fem_obj.fun_analytical(gauss_points);

        % 计算有限元解在高斯点的值（线性时）
        u_fem_gauss = u_fem_left + (u_fem_right - u_fem_left)*(gauss_points - node_left)/h;

        % 计算误差
        error_gauss = abs(u_analytical - u_fem_gauss);%gauss积分中在x1与x2点的f
        % 累加误差
        L1 = L1 + gauss_weights * error_gauss;
        L2 = L2 + gauss_weights * error_gauss.^2;
        
        % 计算Linf误差
        sample_points = linspace(node_left, node_right, sample_num);%单元内采样点

        u_analytical_sample = fem_obj.fun_analytical(sample_points);
        u_fem_sample = u_fem_left + (u_fem_right - u_fem_left)*(sample_points - node_left)/h;
        error_sample = abs(u_analytical_sample - u_fem_sample);
        Linf_element = max(abs(error_sample));%单元内的Linf误差

        if Linf_element > Linf
            Linf = Linf_element; %更新全局Linf误差
        end
    end

    L2=sqrt(L2);%L2范数取平方根

    error=[L1, L2, Linf];%返回L1,L2,Linf误差
end