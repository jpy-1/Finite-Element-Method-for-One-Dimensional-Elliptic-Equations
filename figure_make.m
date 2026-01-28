%figure_make.m
% 绘制有限元解与解析解的对比图像
clear;
clc;
%% 在同一张图上绘制多组 N 的线性有限元解与解析解比较
fem_obj = fem_class();
a = fem_obj.a;
b = fem_obj.b;

N_list = [40,80,160,320];

figure; hold on; box on;
% 画解析解（细网格）
x_dense = linspace(a, b, 2000);
u_exact = fem_obj.fun_analytical(x_dense);
plot(x_dense, u_exact, 'k--', 'LineWidth', 1.5, 'DisplayName', '解析解');

colors = lines(numel(N_list));
for k = 1:numel(N_list)
    N = N_list(k);
    h = (b - a) / N;
    nodes = a:h:b;           % 1 x (N+1)
    u_fem = fun_fem1(N, fem_obj); % 返回 (N+1) x 1 向量
    plot(nodes, u_fem, '-o', 'Color', colors(k,:), 'LineWidth', 1, ...
        'MarkerSize',4, 'DisplayName', sprintf('FEM N=%d', N));
end

xlabel('x'); ylabel('u(x)');
title('分片线性基函数下不同剖分的有限元解与解析解对比');
legend('Location','best');
hold off;


%% 在同一张图上绘制多组 N 的二次有限元解与解析解对比图像
figure; hold on; box on;
% 画解析解（细网格）
plot(x_dense, u_exact, 'k--', 'LineWidth', 1.5, 'DisplayName', '解析解');

colors = lines(numel(N_list));
for k = 1:numel(N_list)
    N = N_list(k);
    h = (b - a) / N;
    nodes=zeros(2*N_list(k)+1,1); %节点
    for i=1:N_list(k)
        nodes(2*i-1)=a+(i-1)*h; %偶数节点
        nodes(2*i)=a+(i-1)*h+h/2; %奇数节点
    end
    nodes(2*N_list(k)+1)=b; %最后一个节点           % 1 x (N+1)
    u_fem = fun_fem2(N, fem_obj); % 返回 (N+1) x 1 向量
    plot(nodes, u_fem, '-o', 'Color', colors(k,:), 'LineWidth', 1, ...
        'MarkerSize',4, 'DisplayName', sprintf('FEM N=%d', N));
end
xlabel('x'); ylabel('u(x)');
title('分片二次多项式基函数下不同剖分的有限元解与解析解对比');
legend('Location','best');
hold off;



%% N=80下线性基函数与二次基函数和解析解的有限元的比较
N_compare = 80; % 选择一个N值进行比较
h = (b - a) / N_compare;
% 线性基函数有限元解
nodes_linear = a:h:b;           % 1 x (N+1)
u_fem_linear = fun_fem1(N_compare, fem_obj); % 返回 (N+1) x 1 向量
% 二次基函数有限元解
nodes_quadratic = zeros(2*N_compare+1,1); %节点
for i=1:N_compare
    nodes_quadratic(2*i-1)=a+(i-1)*h; %偶数节点
    nodes_quadratic(2*i)=a+(i-1)*h+h/2; %奇数节点
end
nodes_quadratic(2*N_compare+1)=b; %最后一个节点           % 1 x (2N+1)
u_fem_quadratic = fun_fem2(N_compare, fem_obj); % 返回 (2N+1) x 1 向量

x_dense = linspace(a, b, 2000);
u_exact = fem_obj.fun_analytical(x_dense);
figure; hold on; box on;
plot(x_dense, u_exact, 'k--', 'LineWidth', 1.5, 'DisplayName', '解析解');
plot(nodes_linear, u_fem_linear, '-o', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1, ...
    'MarkerSize',4, 'DisplayName', '分片线性基函数');
plot(nodes_quadratic, u_fem_quadratic, '-s', 'Color', [0.0000, 0.4470, 0.7410], 'LineWidth', 1, ...
    'MarkerSize',4, 'DisplayName', '分片二次多项式基函数');
xlabel('x'); ylabel('u(x)');
title(sprintf('N=%d 分片线性基函数与分片二次多项式基函数与解析解有限元解对比',N_compare));
legend('Location','best');
hold off;

%% N=160下线性基函数与二次基函数和解析解的有限元的比较
N_compare = 160; % 选择一个N值进行比较
h = (b - a) / N_compare;
% 线性基函数有限元解
nodes_linear = a:h:b;           % 1 x (N+1)
u_fem_linear = fun_fem1(N_compare, fem_obj); % 返回 (N+1) x 1 向量
% 二次基函数有限元解
nodes_quadratic = zeros(2*N_compare+1,1); %节点
for i=1:N_compare
    nodes_quadratic(2*i-1)=a+(i-1)*h; %偶数节点
    nodes_quadratic(2*i)=a+(i-1)*h+h/2; %奇数节点
end
nodes_quadratic(2*N_compare+1)=b; %最后一个节点           % 1 x (2N+1)
u_fem_quadratic = fun_fem2(N_compare, fem_obj); % 返回 (2N+1) x 1 向量
x_dense = linspace(a, b, 2000);
u_exact = fem_obj.fun_analytical(x_dense);
figure; hold on; box on;
plot(x_dense, u_exact, 'k--', 'LineWidth', 1.5, 'DisplayName', '解析解');
plot(nodes_linear, u_fem_linear, '-o', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1, ...
    'MarkerSize',4, 'DisplayName', '分片线性基函数');
plot(nodes_quadratic, u_fem_quadratic, '-s', 'Color', [0.0000, 0.4470, 0.7410], 'LineWidth', 1, ...
    'MarkerSize',4, 'DisplayName', '分片二次多项式基函数');
xlabel('x'); ylabel('u(x)');
title(sprintf('N=%d 分片线性基函数与分片二次多项式基函数与解析解有限元解对比',N_compare));
legend('Location','best');
hold off;


