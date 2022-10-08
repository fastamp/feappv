## 8增强假设应变法
### 8.1.Hu-Washizu变分原理
线弹性问题的变分原理如下：
$$
\begin{aligned}
    \Pi(\boldsymbol{u,\sigma,\epsilon})=&\frac{1}{2}\int_{\Omega}\boldsymbol{\epsilon^TD\epsilon}d\Omega-\int_\Omega\boldsymbol{\epsilon^TD\epsilon^0}d\Omega+\\
    &\int_\Omega\boldsymbol{\sigma^T}(\nabla^{(s)}\boldsymbol{u}-\boldsymbol{\epsilon})d\Omega-\int_\Omega\boldsymbol{u^Tb_v}d\Omega\\
    &-\int_{\Gamma_t}\boldsymbol{u^T\bar{t}}d\Gamma-\int_{\Gamma_u}\boldsymbol{t^T(u-\bar{u})}d\Gamma\\
    =&Stationary\tag{8.1}
\end{aligned}
$$
### 8.4.非线性弹性
非线性超弹性材料的应力直接根据应变能密度函数$W(\boldsymbol{\epsilon})$求导得到：
$$
\boldsymbol{\sigma}=\frac{\partial W}{\partial\boldsymbol{\epsilon}}\tag{8.52}
$$
分量格式：
$$
\sigma_{ij}=\frac{\partial W}{\partial\epsilon_{ij}}\tag{8.53}
$$
对于线弹性材料有：
$$
W(\boldsymbol{\epsilon})=\frac{1}{2}\boldsymbol{\epsilon^TD\epsilon-\epsilon^TD\epsilon^0}\tag{8.54}
$$
对于增强应变公式，应力计算如下：
$$
\tilde{\boldsymbol{\sigma}}=\left.\frac{\partial W}{\partial\boldsymbol{\epsilon}}\right|_{\boldsymbol{\epsilon}=\nabla^{(s)}\boldsymbol{u}+\tilde{\boldsymbol{\epsilon}}}\tag{8.55}
$$
变分方程：
$$
\frac{d\Pi_e}{d\eta}=\int_{\Omega_e}(\nabla^{(s)}\boldsymbol{U})^T\tilde{\boldsymbol{\sigma}}d\Omega+\int_{\Omega_e}\tilde{\boldsymbol{E}}^T\tilde{\boldsymbol{\sigma}}d\Omega\tag{8.56}
$$
残力：
$$
\boldsymbol{R}_I=\boldsymbol{F}_I-\int_{\Omega_e}\boldsymbol{B_I^T\tilde{\sigma}}d\Omega\tag{8.57}
$$
增强模式的残力(在每个单元消失)：
$$
\tilde{\boldsymbol{R}}_\alpha=-\int_{\Omega_e}\boldsymbol{\psi_\alpha^T\tilde{\sigma}}d\Omega=\bold{0}\tag{8.58}
$$
### 8.5求解方法：牛顿法
牛顿法流程：
1.给定方程组
$$
\boldsymbol{f(x)=0}\tag{8.59}
$$
2.构建方程在当前点$\boldsymbol{x^{(i)}}$的线性部分：
$$
\boldsymbol{f^{(i+1)}\simeq f^{(i)}}+\left.\frac{\partial \boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{\boldsymbol{x=x^{(i)}}}d\boldsymbol{x^{(i+1)}=\boldsymbol{0}}\tag{8.60}
$$
其中$d\boldsymbol{x^{(i+1)}}$是$\boldsymbol{x}$的增量。  
3.求解线性问题：
$$
d\boldsymbol{x^{(i+1)}}=-(\boldsymbol{F^{(i)}})^{-1}\boldsymbol{f}^{(i)};\boldsymbol{F^{(i)}}=\left.\frac{\partial \boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{\boldsymbol{x=x^{(i)}}}\tag{8.61}
$$
$\boldsymbol{F^{(i)}}$是方程的雅克比矩阵(切线刚度矩阵)，更新变量：
$$
\boldsymbol{x}^{(i+1)}=\boldsymbol{x}^{(i)}+d\boldsymbol{x}^{(i)}\tag{8.62}
$$
4.重复2和3，直到满足收敛条件：
$$
\left|d\boldsymbol{x}^{(i+1)}\right|\lt tol\left|\boldsymbol{x}^{(i+1)}\right|
$$