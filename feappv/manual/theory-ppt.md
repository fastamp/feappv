## 8��ǿ����Ӧ�䷨
### 8.1.Hu-Washizu���ԭ��
�ߵ�������ı��ԭ�����£�
$$
\begin{aligned}
    \Pi(\boldsymbol{u,\sigma,\epsilon})=&\frac{1}{2}\int_{\Omega}\boldsymbol{\epsilon^TD\epsilon}d\Omega-\int_\Omega\boldsymbol{\epsilon^TD\epsilon^0}d\Omega+\\
    &\int_\Omega\boldsymbol{\sigma^T}(\nabla^{(s)}\boldsymbol{u}-\boldsymbol{\epsilon})d\Omega-\int_\Omega\boldsymbol{u^Tb_v}d\Omega\\
    &-\int_{\Gamma_t}\boldsymbol{u^T\bar{t}}d\Gamma-\int_{\Gamma_u}\boldsymbol{t^T(u-\bar{u})}d\Gamma\\
    =&Stationary\tag{8.1}
\end{aligned}
$$
### 8.4.�����Ե���
�����Գ����Բ��ϵ�Ӧ��ֱ�Ӹ���Ӧ�����ܶȺ���$W(\boldsymbol{\epsilon})$�󵼵õ���
$$
\boldsymbol{\sigma}=\frac{\partial W}{\partial\boldsymbol{\epsilon}}\tag{8.52}
$$
������ʽ��
$$
\sigma_{ij}=\frac{\partial W}{\partial\epsilon_{ij}}\tag{8.53}
$$
�����ߵ��Բ����У�
$$
W(\boldsymbol{\epsilon})=\frac{1}{2}\boldsymbol{\epsilon^TD\epsilon-\epsilon^TD\epsilon^0}\tag{8.54}
$$
������ǿӦ�乫ʽ��Ӧ���������£�
$$
\tilde{\boldsymbol{\sigma}}=\left.\frac{\partial W}{\partial\boldsymbol{\epsilon}}\right|_{\boldsymbol{\epsilon}=\nabla^{(s)}\boldsymbol{u}+\tilde{\boldsymbol{\epsilon}}}\tag{8.55}
$$
��ַ��̣�
$$
\frac{d\Pi_e}{d\eta}=\int_{\Omega_e}(\nabla^{(s)}\boldsymbol{U})^T\tilde{\boldsymbol{\sigma}}d\Omega+\int_{\Omega_e}\tilde{\boldsymbol{E}}^T\tilde{\boldsymbol{\sigma}}d\Omega\tag{8.56}
$$
������
$$
\boldsymbol{R}_I=\boldsymbol{F}_I-\int_{\Omega_e}\boldsymbol{B_I^T\tilde{\sigma}}d\Omega\tag{8.57}
$$
��ǿģʽ�Ĳ���(��ÿ����Ԫ��ʧ)��
$$
\tilde{\boldsymbol{R}}_\alpha=-\int_{\Omega_e}\boldsymbol{\psi_\alpha^T\tilde{\sigma}}d\Omega=\bold{0}\tag{8.58}
$$
### 8.5��ⷽ����ţ�ٷ�
ţ�ٷ����̣�
1.����������
$$
\boldsymbol{f(x)=0}\tag{8.59}
$$
2.���������ڵ�ǰ��$\boldsymbol{x^{(i)}}$�����Բ��֣�
$$
\boldsymbol{f^{(i+1)}\simeq f^{(i)}}+\left.\frac{\partial \boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{\boldsymbol{x=x^{(i)}}}d\boldsymbol{x^{(i+1)}=\boldsymbol{0}}\tag{8.60}
$$
����$d\boldsymbol{x^{(i+1)}}$��$\boldsymbol{x}$��������  
3.����������⣺
$$
d\boldsymbol{x^{(i+1)}}=-(\boldsymbol{F^{(i)}})^{-1}\boldsymbol{f}^{(i)};\boldsymbol{F^{(i)}}=\left.\frac{\partial \boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{\boldsymbol{x=x^{(i)}}}\tag{8.61}
$$
$\boldsymbol{F^{(i)}}$�Ƿ��̵��ſ˱Ⱦ���(���߸նȾ���)�����±�����
$$
\boldsymbol{x}^{(i+1)}=\boldsymbol{x}^{(i)}+d\boldsymbol{x}^{(i)}\tag{8.62}
$$
4.�ظ�2��3��ֱ����������������
$$
\left|d\boldsymbol{x}^{(i+1)}\right|\lt tol\left|\boldsymbol{x}^{(i+1)}\right|
$$