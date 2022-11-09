
# 电子哈密顿量
## 总哈密顿量
分子体系的总哈密顿量表示为
$$
\hat H_{tot} = \hat T_n + \hat T_e + \hat V_{nn} + \hat V_{ne} + \hat V_{ee}
$$
其中包含两类五项：
- 动能：核，电子
- 势能：核-核，电子-电子，核-电子

## BO近似，绝热近似
- 电子速度 >> 核运动速度，核相对静止（核动能为0）
- 核感受的电子势能平均化


## 单电子哈密顿量
$$
\begin{aligned}
\hat H_{el} =& \hat H_{tot} - \hat T_n = \hat T_e + \hat V_{nn} + \hat V_{ne} + \hat V_{ee}
\\=&\sum_i -\frac{1}{2} \nabla_i^2 - \sum_i\sum_A\frac{Z_A}{|r_i-R_A|} + \frac{1}{2}\sum_i\sum_j \frac{1}{|r_i-r_j|} + \hat V_{nn}
\end{aligned}
$$
其中 $-\frac{1}{2} \nabla_i^2$ 为动能算符，$Z_A$ 为第 $A$ 个原子核的电荷量（正数），$R_A$ 为原子核的空间坐标，$r_i$ 为第 $i$ 个电子的空间坐标，第三项为任意两个电子间的排斥能，最后一项为核间互斥能，体系给定它就是一个常数，后续推导可以暂时忽略此项。

将多体问题分解为求解单电子本征值方程:
$$
\hat H_{el} | \Psi_i^{el} \rangle = E_{el} | \Psi_i^{el} \rangle
$$

Hartree-Fock 目的：选择一个最好的组态，使得体系能量最小化。关键在于如何处理平均场即 $\hat V_{ee}$ 。

# 波函数

Slater行列式表示
- 交换反对称
- Pauli不相容原理

双电子实例：
$$
\Psi(\phi_1,\phi_2)=
\begin{vmatrix}
\phi_1(x_1) & \phi_1(x_2)\\
\phi_2(x_1) & \phi_2(x_2)\\
\end{vmatrix} 
= \frac{1}{\sqrt 2}[\phi_1(x_1)\phi_2(x_2)-  \phi_1(x_2)\phi_2(x_1)]
$$
其中 $\phi_i(x_j)$ 表示第 $j$ 个电子在第 $i$ 个分子轨道状态，$x_i=(r_i,w_i)$

# 能量表达式

以2电子体系举例

能量表达式：展开12项，最后剩余4项

一些简化的记号

$$
E = \langle \Psi | \hat H | \Psi \rangle
\\= \sum_i\langle i | h | i\rangle  + \frac{1}{2} \sum_i\sum_j \langle ij||ij\rangle
$$
其中单电子积分记号为：
$$
\langle i | h | j\rangle = \int \phi_i^*(x)h(r)\phi_j(x) dx
$$
双电子积分记号为：
$$
\langle ij|kl\rangle = \iint \phi_i^*(x_1)\phi_j^*(x_2)\frac{1}{|r_1-r_2|}\phi_k(x_1) \phi_l(x_2) dx_1dx_2
\\ \langle ij||ij\rangle = \langle ij|ij\rangle - \langle ij||ji\rangle
$$

# 能量最小化

变分法，用拉格朗日乘子法保证约束条件，分子轨道满足正交归一
限制条件：$\langle i | j \rangle = \delta_{ij} => \langle i | j \rangle - \delta_{ij} = 0$
拉格朗日乘子法：
$$
L = E - \sum_{ij} \lambda_{ij} (\langle i | j \rangle - \delta_{ij})
$$

扰动来获取对应的结果
改变第 m 个轨道，$\phi_m=(1+\delta)\phi_m$，简洁记为 $m=m+\delta m$。
极值条件：$\frac{\partial L}{\partial m}=0$ 等价于

$$
\delta L = L(m=m+\delta m) - L(m) = 0
$$

变分：针对单电子，双电子，拉格朗日系数分别做微分，最后合并结果（仅考虑一阶贡献）。最后得到HF方程：
$$
\hat f_i \phi_m = \sum_i\lambda_{mi} \phi_i
$$

其中 Fock 算符为 $\hat f_i = \hat h_i + \hat J - \hat K$，式子中 $\hat J$ 为库伦算符，$\hat K$ 为交换算符（强制定义，仅仅为了一致性），其具体定义如下：
$$
    \hat J g_m(x_2) = \sum_i \int \phi_i^*(x_1) \frac{1}{r_{12}}\phi_i(x_1) dx_1 g_m(x_2)\\
    \hat K g_m(x_2) = \sum_i \int \phi_i^*(x_1) \frac{1}{r_{12}}\phi_i(x_2) dx_1 g_m(x_1)\\
$$

# 编程实现

## 正则HF方程

通过酉变换和对角化将其变为正则HF方程
选一个能够对角化 $\lambda$ 系数矩阵的酉矩阵 $U$，即 $U\lambda U^{-1}$ 为对角矩阵，从而将波函数转换为 $U\Phi$，其具体推导如下：
$$
F\Phi = \lambda \Phi\\
UF\Phi = U\lambda U^{-1}U \Phi\\
F(U\Phi) = (U\lambda U^{-1})(U \Phi)\\
F\Phi' = \epsilon \Phi'\\
$$
因此能够得到正则HF方程：
$$
\hat f_i \phi_i = \epsilon_i \phi_i
$$
标准的本征值问题。

## Roothaan-Fall 方程

若想在计算机上实现，必须将分子轨道用原子轨道展开，即分子轨道是原子轨道的线性组合：$\phi_i=\sum_uC_{ui}\chi_u $。将该式子代入正则HF方程，并同时左乘一个基底函数并积分，即可得到Roothaan-Fall 方程：

$$
FC=SCE
$$
其中 
- F是 Fock 矩阵：$F_{uv} = \langle \chi_u| \hat f | \chi_v \rangle=h_{uv}+\sum_{\alpha\beta}D_{\alpha\beta}[\langle \alpha u|\beta v\rangle -\langle \alpha u|v\beta\rangle]$，其中密度矩阵 $D_{uv}=\sum_iC_{ui}^*C_{vi}$
- S是重叠矩阵：$S_{uv}=\langle \chi_u | \chi_v \rangle$，原子轨道不保证正交归一
- C是系数矩阵
- E是对角矩阵，代表每个分子轨道的能量

## 自洽场求解（SCF）
存在循环依赖：F->D->C->F

迭代求解算法

- 给定原子基组 $\{\chi_i\}$，计算重叠矩阵 $S$，$h_{uv}$，双电子积分
- 猜一个 $C^{(0)}$，进而计算出 $D^{(0)},F^{(0)}$
- 对角化 $FC=SCE$，得到 $C^{(1)}$
- 重复第2,3步骤，直至收敛（D前后两次变量低于某个阈值）

对角化方程转换

# References
对于量子化学初学者（原计算机选手，数学基础仅有高数，线代，概率论，离散数学）想要快速入门理解，我采用的路线是忽略大量的基础数学知识（一开始就看这个很容易觉得无聊从而放弃，别问我咋知道的），直接看hartree-fock的推导视频，在看视频中将不懂得基础概念和数学方法粗略了解下（先掌握脉络，后面再深入），然后跟着推公式，自己整理笔记，发现不通畅时回去看视频或者看书/文章，解决一个个的问题。最后是将公式转换为代码，来检验自己的推导是否正确，是否真的理解学习内容。

这么一套下来，在学中练，练中学，理论与实践并进，同时目标明确，反馈及时，能够促进思考与实践效率。

## 理论推导
[量子化学入门【已完结】](https://www.bilibili.com/video/BV12k4y1d73q/?p=21&share_source=copy_web&vd_source=8f1a8bf0da41ecd44f5ab605e6f4e43b)：白板推导教学，强烈推荐新手小白观看学习

[Hartree-Fock超详细推导](https://www.cnblogs.com/miccoui/p/15410545.html)：文字版，非常详细，可对照视频看，但不包括自旋部分

[从Hartree-Fock方程推导Roothaan方程](https://www.bilibili.com/video/BV1u14y157Af/?share_source=copy_web&vd_source=8f1a8bf0da41ecd44f5ab605e6f4e43b)：白板推导，包括了自旋部分，能得到严格的RHF的求解公式

《Modern quantum chemistry introduction to advanced electronic structure theory》 Attila Szabo，适合入门观看，HF部分教你实现Hartree-Fock代码

《量子化学中册》徐光宪，非常详实，想了解数学和技术细节的冲他，但初学最好和Szabo，上面推荐的视频一块食用

《Quantum Chemistry》：Levine，比较适合入门选手

## 代码实践

[hf-tutorial](https://github.com/yangdatou/hf-tutorial)：提供一个测试框架，让你填空实现RHF，各类积分，重叠矩阵，哈密顿矩阵都作为输入，你无需关心，只要实现最核心的SCF算法，同时提供了一份已有的解决方案作为参考，非常适合检验自己的学习成果。

[Psi4Numpy](https://github.com/psi4/psi4numpy)：基于Psi4和Numpy实现各类量化程序的教程

[xDH在python下实现的简易教程](https://py-xdh.readthedocs.io/zh_CN/latest/index.html)：量化介绍及其实现原理实践