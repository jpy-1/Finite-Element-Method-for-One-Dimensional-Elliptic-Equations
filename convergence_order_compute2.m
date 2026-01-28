% convergence_order_compute2.m 计算使用二次基函数空间的有限元解的收敛阶
clear; clc;
fem_obj=fem_class();
N_list=[20,40,80,160,320];
N_initial=N_list(1);%初始网格数
num_N=length(N_list);
% 初始化存储结果的矩阵
results = zeros(num_N, 7);  % 列依次为：N, L1误差, L1阶,L2误差, L2阶,Linf误差, Linf阶
results(:,1)=N_list';

u_fem_1=fun_fem2(N_initial,fem_obj);
error_1=error_compute2(u_fem_1,fem_obj,N_initial);
L1_error_1=error_1(1);
L2_error_1=error_1(2);
Linf_error_1=error_1(3);
results(1,2:2:6)=[L1_error_1, L2_error_1, Linf_error_1];
disp('二次基函数的有限元解随着网格加密误差变化情况：');
fprintf('当N=%d, L1误差: %e, L2误差: %e, Linf误差: %e\n', N_initial, L1_error_1, L2_error_1, Linf_error_1);


for i=2:num_N
    N=N_list(i);
    %误差计算
    u_fem_2=fun_fem2(N,fem_obj);
    error_2=error_compute2(u_fem_2,fem_obj,N);
    L1_error_2=error_2(1);
    L2_error_2=error_2(2);
    Linf_error_2=error_2(3);
    fprintf('当N=%d, L1误差: %e, L2误差: %e, Linf误差: %e\n', N, L1_error_2, L2_error_2, Linf_error_2);

    %收敛阶的计算
    L1_order=log(L1_error_1/L1_error_2)/log(2);
    L2_order=log(L2_error_1/L2_error_2)/log(2);
    Linf_order=log(Linf_error_1/Linf_error_2)/log(2);
    fprintf('当N=%d, L1收敛阶: %e, L2收敛阶: %e, Linf收敛阶: %e\n', N, L1_order, L2_order, Linf_order);

    %结果存储
    results(i,2:end)=[L1_error_2, L1_order, L2_error_2, L2_order, Linf_error_2, Linf_order];


    L1_error_1=L1_error_2;
    L2_error_1=L2_error_2;
    Linf_error_1=Linf_error_2;
end

% disp('-----------------------------------------------------');
% disp('网格数 N |    L1误差   | L1阶  |    L2误差   | L2阶  |   Linf误差  | Linf阶  |');
% disp('-----------------------------------------------------');
% disp(results);

% 将结果写入 CSV 文件（数值按 6 位有效数字）
filename = fullfile(pwd, 'error_convergence2.csv'); % 保存到当前工作目录

% 准备带表头的字符串表格，数值按 6 位有效数字格式化
headers = {'N','L1_error','L1_order','L2_error','L2_order','Linf_error','Linf_order'};
num_rows = size(results,1);
formatted = cell(num_rows, numel(headers));
for i = 1:num_rows
    formatted{i,1} = sprintf('%d', results(i,1)); % N 用整数格式
    for j = 2:size(results,2)
        formatted{i,j} = sprintf('%.7g', results(i,j)); % 7 位有效数字
    end
end

T = cell2table(formatted, 'VariableNames', headers);
writetable(T, filename);
fprintf('结果已写入 %s\n', filename);

