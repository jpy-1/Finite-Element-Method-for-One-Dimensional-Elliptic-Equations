# Finite-Element-Method-for-One-Dimensional-Elliptic-Equations
## Using polynomial functions as basis functions, including both linear and quadratic element methods, and taking the one-dimensional elliptic equation as an example, the error calculation and convergence order verification are computed and presented.
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
