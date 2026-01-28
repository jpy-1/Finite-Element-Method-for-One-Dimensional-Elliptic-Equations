% error_compute2.m 计算使用二次基函数空间的有限元解的L1, L2, Linf误差
function error = error_compute2(u_fem,fem_obj,N)
    % 初始化误差
    L1 = 0;
    L2 = 0;
    Linf = 0;
    a=fem_obj.a; %左端点
    b=fem_obj.b; %右端点
    h=(b-a)/N; %步长
    nodes=zeros(2*N+1,1); %节点
    for i=1:N
        nodes(2*i-1)=a+(i-1)*h; %偶数节点
        nodes(2*i)=a+(i-1)*h+h/2; %奇数节点
    end
    nodes(2*N+1)=b; %最后一个节点
    sample_num=50;%L^\infty误差单元内采样点数量

    %遍历每个单元计算误差
    for i=1:N
        % 当前单元e_i: [x_{i-1}, x_i]
        node_left = nodes(2*i-1);
        node_mid = nodes(2*i);
        node_right = nodes(2*i+1);
        u_fem_left = u_fem(2*i-1);       % 单元左端点的有限元解
        u_fem_mid = u_fem(2*i);         % 单元中点的有限元解
        u_fem_right = u_fem(2*i+1);    % 单元右端点的有限元解

        %使用Gauss积分计算L1,L2误差
        gauss_points = [node_left + (1-1/sqrt(3))*h/2; node_left + (1+1/sqrt(3))*h/2];
        gauss_weights = [h/2, h/2];

        % 计算解析解在高斯点的值
        u_analytical = fem_obj.fun_analytical(gauss_points);

        % 计算有限元解在高斯点的值（二次时）
        gauss_matrix=[(sqrt(3)+1)/6,(1-sqrt(3))/6;
                        2/3, 2/3;
                        (1-sqrt(3))/6,(1+sqrt(3))/6];% 二次基函数在高斯点的值
        u_fem_gauss = gauss_matrix'*[u_fem_left;u_fem_mid; u_fem_right];

        %计算误差
        error_gauss = abs(u_analytical - u_fem_gauss);%gauss积分中在x1与x2点的f
        % 累加误差
        L1 = L1 + gauss_weights * error_gauss;
        L2 = L2 + gauss_weights * error_gauss.^2;
        
        % 计算L^\infty误差，单元内均匀采样
        sample_points = linspace(node_left, node_right, sample_num);
        u_analytical_samples = fem_obj.fun_analytical(sample_points);
        u_fem_samples= u_fem_left * ((sample_points - node_mid).*(sample_points - node_right))/(h^2/2) + ...
                        u_fem_mid * ((sample_points - node_left).*(sample_points - node_right))/(-h^2/4) + ...
                        u_fem_right * ((sample_points - node_left).*(sample_points - node_mid))/(h^2/2);
        error_samples = abs(u_analytical_samples - u_fem_samples);
        Linf_element = max(abs(error_samples));%单元内的Linf误差

        if Linf_element > Linf
            Linf = Linf_element; %更新全局Linf误差
        end
    end

    L2=sqrt(L2);%L2范数取平方根
    
    error=[L1, L2, Linf];%返回L1,L2,Linf误差
end
