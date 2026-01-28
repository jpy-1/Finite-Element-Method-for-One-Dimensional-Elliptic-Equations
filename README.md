# Finite-Element-Method-for-One-Dimensional-Elliptic-Equations
## 采用多项式函数作为基函数，涵盖线性元与二次元两种方法，并以一维椭圆方程为例，完成了误差计算与收敛阶验证，并给出相应结果。
## The following is the solution to the equation：
### equation
$$
-u'' = \sin 4x,\quad x\in \left[0,\frac{\pi}{8}\right]
$$

### boundary condition
$$
u(0) = 0,\quad u'\left(\frac{\pi}{8}\right) = 0
$$

### finite space
$$
S\big|_{[x_{i-1},x_i]} \in P^k
$$
## 代码使用：后缀1,2分别表示一次元与二次元
### fem_class为求解方程，可在其中修改参数更改方程
### fun_fem为求解函数，输入求解方程，输出数值解
### errror_compute函数，输入数值解与解析解，给出有限元空间的 $$L^1,L^2,L^{\infty}$$ 误差
### convergence_order_compute函数：该程序的主函数，调用其余函数给出误差，并计算不同误差下的收敛阶
### Numerical Experiment Report.pdf 数值实验报告，包含求解图像与误差表格
