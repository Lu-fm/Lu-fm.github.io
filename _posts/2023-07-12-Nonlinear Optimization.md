---
layout:     post
title:      "非线性优化概述"
# subtitle:   " "\The base of everything flying\""
date:       2023-07-12 19:51:00
author:     "lufangmin"
header-img: "img/shuttle.jpg"
catalog: true
tags:
    - Optimization
    - Direct Method
---
uploaded on 01/09/24

### 1 线搜索法
#### 1.1 简介
线搜索法的迭代过程如下：

$$
x_{k+1} = x_{k}+\alpha_{k}p_{k}
$$

其中$p_{k}$是搜索方向，$\alpha_{k}$是步长。

大多数线搜索法要求$p_{k}$是下降方向，即：

$$
p_{k}^{T}\nabla f_{k}<0
$$

搜索方向通常有如下形式：

$$
p_{k} = -B_{k}^{-1}\nabla f_{k}
$$

其中$B_{k}$是非奇异的对称阵。当$B_{k}$是正定矩阵时，有

$$
p_{k}^{T}\nabla f_{k} = -p_{k}^{T}B_{k}^{-1}\nabla f_{k}<0
$$

因此$p_{k}$是下降方向。

当$B_{k}$是单位阵时，这是最速下降法；
当$B_{k}$是Hessian矩阵（$\nabla^{2}f(x_{k})$）时，这是牛顿法；
当$B_{k}$是近似Hessian矩阵时，这是拟牛顿法。

#### 1.2 步长选择

步长选择过小，可能满足下降的要求，但是无法迭代到局部最小值，因此我们需要**充分下降**。

1）Armijo condition
$$
f(x_{k}+\alpha_{k}p_{k})\leq f(x_{k})+c_{1}\alpha_{k}\nabla f_{k}^{T}p_{k}
$$
其中$c\in(0,1)$。

![[Pasted image 20230712103740.png]]

2）Wolfe conditions（防止步长过小）
$$
\begin{align}
f(x_{k}+\alpha_{k}p_{k})&\leq f(x_{k})+c_{1}\alpha_{k}\nabla f_{k}^{T}p_{k} \\
\nabla f(x_{k}+\alpha_{k}p_{k})^{T}p_{k}&\geq c_{2}\nabla f_{k}^{T}p_{k}\
\end{align}
$$
其中$0<c_{1}<c_{2}<1$。
![[Pasted image 20230712104424.png]]

对任何光滑且有下界的$f$， 可证明存在$\alpha$满足Wolfe conditions.

3) Goldstein conditions
$$
\begin{align}
f(x_{k})+(1-c)\alpha_{k}\nabla f_{k}^{T}p_{k}&\leq f(x_{k}+\alpha_{k}p_{k}) \\
&\leq
f(x_{k})+c\alpha_{k}\nabla f_{k}^{T}p_{k} 
\end{align}
$$
其中$0<c<1/2$
![[Pasted image 20230712105013.png]]


> **Algorigthm**  (Backtracking Line Search)
> choose $\overline{\alpha}>0,\rho\in(0,1), c\in(0,1); \text{set } \alpha \leftarrow \overline{\alpha}$
> 
> repeat until $f(x_{k}+\alpha_{k}p_{k})\leq f(x_{k})+c_{1}\alpha_{k}\nabla f_{k}^{T}p_{k}$
> 		$\alpha\leftarrow \rho \alpha$
> end
> Terminate with $\alpha_{k}=\alpha$

where $\rho\in (0,1)$

### 1.3 步长选择方法

1）$f(x)$ 是凸的二次型， 即$f(x) = x^{T}Qx-b^{T}x$, 则
$$
\alpha_{k} = \mathop{\arg \min}_{\alpha}f(x_{k}+\alpha p_{k})=-\frac{\nabla f_{k}^{T}P_{k}}{p_{k}^{T}Qp_{k}}
$$

2）插值，使用二次或三次多项式（代理函数）进行插值，使得其满足：
$$
\begin{align}
\phi_{c}(0)&=\phi(0), \phi_{c}'(0)=\phi'(0),  \\
\phi_{c}(\alpha_{k})&=\phi(\alpha_{k}), \phi_{c}(\alpha_{k-1})=\phi(\alpha_{k-1})
\end{align}

$$
令$\phi_{c}$导数为0即可得到$\alpha_{k+1}$


### 2 信赖域法

线搜索法基于二次函数模型得到搜索方向，再寻找合适的步长；信赖域法同时确定搜索方向以及相应的步长。

我们介绍以如下最优化子问题为基础的信赖域法

$$
\begin{align}
&\min_{p\in R^{n}} m_{k}(p)=f_{k}+g_{k}^{T}p+\frac{1}{2}p^{T}B_{k} \\
&s.t. \quad \Vert p\Vert \leq \Delta_{k}

\end{align}
$$

其中$B_{k}$是对称阵，当其为Hessian$\nabla^{2}f_{k}$时，该方法为牛顿信赖域法。$\Delta_{k}$是信赖域半径。

利用拉格朗日乘子法可得到该子问题的解的KKT条件：
$$
\begin{align} 
(B+\lambda I)p^{\ast}&=-g \\
\lambda(\Delta- \Vert p^{\ast}\Vert)&=0 \\
(B+\lambda I)&\geq0
\end{align}
$$



#### 2.1 Cauchy Point Calculation

1）首先忽略二次项，求解如下线性子问题
$$
p_{k}^{s}= \mathop{\arg\min}_{p\in R^{n}} \quad f_{k}+g_{k}^{T}p\qquad

s.t. \quad \Vert p\Vert \leq \Delta_{k}

$$
得到
$$
p_{k}^{s} = -\frac{\Delta_{k}}{||g_{k}||}g_{k}
$$
再求解标量$\tau$满足：
$$
\tau_{k}= \mathop{\arg\min}_{\tau\geq_{0}} \quad m_{k}(\tau_{k}p_{k}^{s}) \qquad

s.t. \quad \Vert \tau p_{k}^{s}\Vert \leq \Delta_{k}

$$
易得：
$$
\begin{equation}
\tau_{k} = \begin{cases}
1 & \text{if }g_{k}^{T}Bg_{k}\leq_{0} \\
\min \left( \frac{||g_{k}||^{3}}{\Delta_{k}g_{k}^{T}Bg_{k}}, 1 \right),  & \text{otherwise}
\end{cases}
\end{equation}

$$

最后得到解$p_{k}^{c}= \tau_{k} p_{k}^{s}$， 实际上这种方法和最速下降没太大区别。

#### 2.2 THE DOGLEG METHOD

![[Pasted image 20230712150234.png]]

该方法将最速下降方向($-\nabla f$)和牛顿法($-B^{-1}\nabla f$)相结合， 该方法通常在$B$为正定时适用。

令$p^{U} = -\frac{g^{T}g}{g^{T}Bg}g$, $p^{B} = -B^{-1}g$,  可将搜索路径写为

$$
\begin{equation}
p(\tau) = \begin{cases}
\tau p^{U},  &   \quad 0\leq\tau \leq 1 \\
p^{U} + (\tau-1)(p^{B}-p^{U}),  & \quad1\leq \tau\leq 2
\end{cases}
\end{equation}
$$
易证明$||p(\tau)||$随$\tau$递增， $m(p(\tau))$ 随$\tau$递减，因此，若$||p^{B}||\leq \Delta$, 则$p_{k}=p^{B}$，否则$\tau$由下式决定
$$
|| p^{U}+(\tau-1)(P^{B}-P^{U})||^{2} = \Delta^{2}
$$

#### 2.3 TWO-DIMENSIONAL SUBSPACE MINIMIZATION

该方法将一维子空间的dogleg方法拓展成二维子空间，即令
$$
p\in \text{span} [g, B^{-1}g], \quad||p||\leq \Delta
$$

当$B$不是正定时，令
$$
p\in \text{span} [g, (B+\alpha I)^{-1}g], \quad \alpha\in(-\lambda_{1}, -2\lambda_{1})
$$
其中$\lambda_{1}$是B的负特征值中最小的一个。

当$||(B+\alpha I)^{-1}g\leq \Delta$时，直接取
$$
p = -(B+\alpha I)^{-1}g+v
$$
其中$v^{T}(B+\alpha I)v\leq0$


### 3 共轭梯度法

该方法是求解大规模线性方程组的最有效的方法之一，同时经过修正也能够用于求解非线性优化问题。

### 3. 1 线性共轭梯度法

待求解问题为线性方程组
$$
Ax=b
$$
其中$A$是$n\times n$的对称正定矩阵，该问题也可等价为
$$
\min \phi(x) = \frac{1}{2}x^{T}Ax - b^{T}x
$$
定义残差 $r(x)$
$$r(x)  = \nabla \phi(x) = Ax-b$$

非0向量$\{p_{0}, p_{1}, \cdots, p_{n-1}\}$ 被称为是关于$A$共轭的，其中 $A$ 是对称且正定的矩阵。
$$
p_{i}^{T}Ap_{j}=0, \quad \text{for all }i\neq j 
$$
易证满足这样性质的向量是线性无关的。按如下策略迭代可以证明， 至多迭代$n$步，即可收敛到$\phi(\cdot)$的最优解$x^{\ast}$。

$$
x_{k+1}=  x_{k}+\alpha_{k}p_{k}
$$
$$
\alpha_{k} = -\frac{r_{k}p_{k}}{p_{k}^{T}Ap_{k}}
$$
**证明：**

由于$\{p_{0}, p_{1}, \cdots, p_{n-1}\}$是线性无关的，故其可张成$n$维线性子空间。若$x^{\ast}$为最优解，$x_{0}$为起始点， 则有
$$
x^{\ast}-x_{0} = \sigma_{0}p_{0}+\sigma_{1}p_{1}+\cdots+\sigma_{n-1}p_{n-1}
$$

两边左乘 $p_{k}^{T}A$, 可得
$$
\sigma_{k} = \frac{p_{k}^{T}A(x^{\ast}-x_{0})}{p_{k}^{T}Ap_{k}}
$$

由迭代策略可知， 
$$
x_{k} = x_{0} + \alpha_{0}p_{0}+\cdots+\alpha_{k-1}p_{k-1}
$$
两边同乘$p_{k}^{T}A$, 可得
$$
p_{k}^{T}A(x_{k}-x_{0}) = 0
$$

故有
$$
p_{k}^{T}A(x^{\ast}-x_{0})  = p_{k}^{T}A(x^{\ast}-x_{k})=p_{k}^{T}(b-Ax_{k}) = -p_{k}^{T}r_{k}
$$
回代可知$\sigma_{k} = \alpha_{k}$, 证毕。


正交性质：
$$
r_{k}^{T}p_{i} = 0,\quad i = 0, 1, \cdots,k-1
$$
$$
r_{k+1} = r_{k} + \alpha_{k}Ap_{k}
$$

### 4 拟牛顿法

- SR1
- SR2(BFGS)


### 5 大规模无约束优化



### 5X 拉格朗日对偶问题

设原问题(primal)为:
$$
\begin{align}
&\min f(x) \\
&s.t. c(x)\leq 0 
\end{align}
$$
拉格朗日函数为：
$$
L(x, \lambda) = f(x) - \lambda c(x)
$$
其对偶性能指标为：
$$
q(\lambda) = \inf_{x} L(x,\lambda)
$$
由于$q(\lambda)$取值可能会无穷小， 因此定义域为：
$$
\mathcal{D}=\{\lambda|q(\lambda)>-\infty\}
$$
对偶问题为：
$$
\max_{\lambda} q(\lambda), \quad s.t.\lambda\geq0
$$
易证明$q(\lambda)$是凹函数， $\mathcal{D}$ 是凸集。

1） Weak Duality

对任意$x, \lambda$在各自的可行域内，有:
$$
q(\lambda)\leq f(x)
$$

### 5Y 线性规划

单纯形法（Simplex Method)

2)内点法（Primal-Dual Interior Point Method)

考虑问题
$$
\begin{align}
&\min_{x} c^{T}x \\
& s.t. Ax = b, x\geq 0
\end{align}
$$
拉格朗日函数为
$$
L = c^{T}x - (Ax-b)^{T}\lambda - x^{T}s
$$
KKT条件为
$$
\begin{align}
A^{T}\lambda+s &= c \\
Ax &= b \\
x_{i}s_{i} &= 0 \\
(x,s)&>0
\end{align}
$$

求解满足KKT条件的点，即为求解如下问题：
$$
F(x,\lambda,s) = \begin{bmatrix}
A^{T}\lambda + s - c \\
Ax-b \\
XS1
\end{bmatrix} = 0
$$
其中$X = diag(x_{1}, \cdots, x_{n}),  S = diag(s_{1},\cdots, s_{n})$。

使用牛顿法求解，即为

$$
J_{F}d = -F
$$
具体为如下表达式：
$$
\begin{bmatrix}
0 & A^{T} & I \\
A & 0 & 0 \\
S & 0 & X 
\end{bmatrix}\begin{bmatrix}
\Delta x \\
\Delta \lambda \\
\Delta s
\end{bmatrix} = -\begin{bmatrix}
A^{T}\lambda + s - c \\
Ax-b \\
XS1
\end{bmatrix} 
$$

此为单纯的牛顿法，为满足$(x,s)\geq0$, 步长可能要取得很小，因此算法效率较低。我们将原问题修改为以下问题：

$$
\begin{align}
&\min_{x} c^{T}x - \tau\Sigma_{i=1}^{n}\ln x_{i} \\
& s.t. Ax = b
\end{align}
$$

该问题是严格凸的，所以极小值唯一，且满足LICQ条件，因此拉格朗日函数的解也唯一。

KKT条件变为
$$
F(x,\lambda,s) = \begin{bmatrix}
A^{T}\lambda + s - c \\
Ax-b \\
XS1
\end{bmatrix} = \begin{bmatrix}
0 \\
0 \\
\tau 1
\end{bmatrix}
$$
其中$\tau\geq 0$ ，则牛顿法求解搜索方向的式子为
$$
\begin{bmatrix}
0 & A^{T} & I \\
A & 0 & 0 \\
S & 0 & X 
\end{bmatrix}\begin{bmatrix}
\Delta x \\
\Delta \lambda \\
\Delta s
\end{bmatrix} = -\begin{bmatrix}
A^{T}\lambda + s - c \\
Ax-b \\
XS1-\tau 1
\end{bmatrix} 
$$

在实际算法中，我们令$\tau = \sigma \mu$， 其中$\sigma\in[0,1]$是人为指定的参数，$\mu = \frac{1}{n}x^{T}s$ 是“duality measure“，为对偶问题的解与原问题的之差的均值。


### 6 二次规划

设二次规划问题为
$$
\begin{align}
&\min Q(x) = \frac{1}{2}x^{T}Hx+c^{T}x \\
s.t.\quad &a_{i}^{T}x = b_{i},\quad i=1,\cdots,m_{e} \\
&a_{i}^{T}x\geq b_{i}, \quad i=m_{e}+1,\cdots,m
\end{align}
$$

>引理
>设$x^{\ast}$是上述二次规划问题的局部极小点，则$x^{\ast}$也必是问题
$$
\begin{align}
&\min Q(x) = \frac{1}{2}x^{T}Hx+c^{T}x \\
s.t.\quad &a_{i}^{T}x = b_{i},\quad i\in \mathcal{A}(x^{\ast})
\end{align}
$$
>的局部极小点。反之，若$x^{\ast}$是原问题的可行点，且是上述问题的KT点，并且相应的拉格朗日乘子满足
$$
\lambda^{\ast}\geq 0, \quad i\in \mathcal{I}(x^{\ast})
$$
>则$x^{\ast}$也是原问题的KT点。

其中$\mathcal{A}(x) = \{i\in \mathcal{E}\cup\mathcal{I}|a_{i}^{T}x = b_{i}\}$为active set。

1）积极集法（Active Set Method)

积极集法是一个可行点方法，即每个迭代点都要求是可行点。每次迭代求解一个等式约束的二次规划。

等式约束子问题SQ：
$$
\begin{align}
&\min_{d\in \mathbb{R}^{n}} c^{T}(x_{k}+d)+\frac{1}{2}(x_{k}+d)^{T}H(x_{k}+d) \\
&s.t. a_{i}^{T}d = 0, i\in A(x_{k})
\end{align}
$$

>Algorithm 
>1. 给出可行点$x_{1}$, 令$W_{1}=A(x_{1}), k=1$
>2. 求解子问题SQ得到$d_{k}$, 
>若$d_{k}=0$, 
>	根据式 $\sum_{i\in W_{k}}a_{i}\lambda_{i}=Hx_{k}+c$ 计算$\lambda_{i}$ 
>	若$\lambda_{i}\geq0(i\in W\cap \mathcal{I})$, 则停止迭代；
>	否则：$j=\mathop{\arg\min_{j\in W\cap \mathcal{I}} \lambda_{j}}; x_{k+1}=x_{k};W_{k+1}=W_{k}\backslash \{j\}$
>	
>否则($d_{k}\neq 0)$：
>	根据下式计算$\alpha$
$$
\alpha_{k} = \min \left( 1, \min_{i\notin W_{k}, a_{i}^{T}d_{k}\leq 0} \frac{b_{i}-a_{i}^{T}x_{k}}{a_{i}^{T}d_{k}} \right)
$$
>	$x_{k+1} = x_{k}+\alpha_{k}d_{k}$
>	若有约束不满足(即$x_{k+1}$为不可行点)，则将该约束添加到$W_{k}$中，得到$W_{k+1}$。
>	否则$W_{k+1} = W_{k}$


### 7 罚函数法和增广拉格朗日法

### 7.1 罚函数法
原问题为

$$
\min_{x} f(x) \quad s.t. c_{i}(x) = 0, i\in \mathcal{E},\quad c_{j}(x)\geq0, j\in\mathcal{I}
$$

构造二次罚函数
$$
Q(x;\mu) = f(x) + \frac{\mu}{2}\sum_{i\in\mathcal{E}}c_{i}^{2}(x)+\frac{\mu}{2}\sum_{i\in\mathcal{I}}([c_{i}(x)]^{-1})^{2}
$$

其中$\mu>0$, 是惩罚参数，$[y]^{-1}=\max(-y,0)$。当$\mu\to\infty$时，无约束问题$Q(x,\mu)$的极小值点即为原等式约束问题的极小值点。当问题中含不等式约束时，即$\mathcal{I}\neq \emptyset$，$Q(x,\mu)$可能没有连续的二阶导数。

### 7.2 增广拉格朗日函数法

$$
L_{A}(x,\lambda;\mu) = f(x) - \sum_{i\in\mathcal{E}}\lambda_{i}c_{i}(x) + \frac{\mu}{2}\sum_{i\in\mathcal{E}}c_{i}^{2}(x)
$$

一阶必要条件为
$$
0\approx \nabla L_{A}(x_{k},\lambda_{k};\mu_{k}) = \nabla f(x_{k}) - \sum_{i\in\mathcal{E}}[\lambda_{i}^{k}-\mu_{k}c_{i}(x_{k})]\nabla c_{i}(x_{k})
$$
与原问题对比可得
$$
\lambda^{\ast}_{i} \approx \lambda_{i}^{k}-\mu_{k}c_{i}(x_{k})
$$
这启发我们可以按照下式更新$\lambda_{k}$:
$$
\lambda_{i}^{k+1} = \lambda_{i}^{k} - \mu_{k}c_{i}(x_{k})
$$
算法流程即为在每次迭代中固定$\mu_{k}$, 求解$x_{k} = \mathop{\arg\min} L_{A}(x,\lambda^{k};\mu^{k})$, 更新$\lambda_{i}^{k+1}$, 增大$\mu_{k+1}$, 继续下次迭代。

算法实现：（若存在不等式约束，则引入松弛变量将其转化为等式约束，松弛变量添加边界约束，即$l<s<u$）


### 8 序列二次规划法

#### 8.1 算法介绍
这里介绍的序列二次规划法属于积极集法，可以在线搜索法或信赖域法的框架下实现。

1)考虑一般的非线性等式约束问题：
$$
\min f(x),\quad s.t.\space c(x)=0
$$

记
$$
A(x)^{T} = [\nabla c_{1}(x), \nabla c_{2}(x),\cdots,\nabla c_{m}(x)]
$$

该问题的KKT条件为
$$
F(x,\lambda) = \begin{bmatrix}
\nabla f(x) - A(x)^{T}\lambda \\
c(x)
\end{bmatrix} = 0
$$
$F(x,\lambda)$的雅可比矩阵为
$$
F'(x,\lambda) = \begin{bmatrix}
\nabla^{2}L(x,\lambda) & -A(x)^{T} \\
A(x) & 0 
\end{bmatrix}
$$

牛顿法迭代步长求解：
$$
\begin{bmatrix}
\nabla^{2}L(x,\lambda) & -A(x)^{T} \\
A(x) & 0 
\end{bmatrix}\begin{bmatrix}
p_{k} \\
p_{\lambda}
\end{bmatrix} = -\begin{bmatrix}
\nabla f(x) - A(x)^{T}\lambda \\
c(x)
\end{bmatrix} 
$$
其中，
$$
\begin{bmatrix}
p_{k} \\
p_{\lambda} 
\end{bmatrix}= \begin{bmatrix}
x_{k+1} \\
\lambda_{k+1}
\end{bmatrix}-\begin{bmatrix}
x_{k} \\
\lambda_{k}
\end{bmatrix}
$$

>Assumptions
>(a) A(x) 行满秩
>(b) Hessian矩阵$\nabla^{2}L(x,\lambda)$在约束的切空间是正定的，即
>$d^{T}\nabla^{2}L(x,\lambda)d>0$, 对任意$d\neq 0, A(x)d=0$。

(a)是约束的线性独立调条件（LICQ），(b)条件成立的情况是在极小值点满足二阶充分条件，且$(x,\lambda)$充分靠近极小值点。若假设(a)(b)均满足，可知该迭代能以二阶速度收敛。

牛顿迭代法求解步长也可以看作是求解如下问题的KKT条件：
$$
\begin{align}
&\min_{p} f_{k}+\nabla f^{T}_{k}p + \frac{1}{2}p^{T}\nabla^{2}L_{k}p \\
&s.t. \quad A_{k}p+c_{k} = 0
\end{align}
$$
设该问题的解及拉格朗日乘子为$(p_{k}, l_{k})$, 则KKT条件为
$$
\begin{bmatrix}
\nabla^{2}L(x,\lambda) & -A(x)^{T} \\
A(x) & 0 
\end{bmatrix}\begin{bmatrix}
p_{k} \\
l_{k}
\end{bmatrix} = -\begin{bmatrix}
\nabla f(x) \\
c(x)
\end{bmatrix} 
$$
与牛顿法比较可得$l_{k} = \lambda_{k+1}$。
实际上，根据约束$A_{k}p+c_{k}=0$ 可知 $\nabla_{x} L^{T}p=\nabla f_{k}^{T}p+\lambda_{k}^{T}c_{k}$, 由于$\lambda_{k}^{T}c_{k}$是常数，因此优化问题可写为
$$
\begin{align}
&\min_{p} L_{k}+\nabla_{x} L^{T}_{k}p + \frac{1}{2}p^{T}\nabla^{2}_{xx}L_{k}p \\
&s.t. \quad A_{k}p+c_{k} = 0
\end{align}
$$
即为$L(x_{k},\lambda_{k})$的二次展开近似，约束为线性近似。

2)考虑不等式约束问题

$$
\begin{align}
&\min f(x) \\
 s.t.\space &c_{i}(x)=0, \quad i\in \mathcal{E}  \\
&c_{i}(x) \geq 0, \quad i\in \mathcal{I}
\end{align}
$$

在迭代求解时我们求解如下子问题：
$$
\begin{align}
\min_{p}\space f_{k}+\nabla f_{k}^{T}+&\frac{1}{2}p^{T}\nabla^{2}_{xx}L_{k}p \\
s.t.\quad \nabla c_{i}(x_{k})^{T}p+c_{i}x_{k}&=0,\quad i\in \mathcal{E} \\
\quad \nabla c_{i}(x_{k})^{T}p+c_{i}x_{k}&\geq0,\quad i\in \mathcal{I} 
\end{align}
$$
即可用二次规划的方法去求解子问题。

上述方法称为IQP(Inequality-constrained QP), 除此之外，还有EQP(Equality-constrained QP)


#### 8.2 问题及解决方法

1）elastic mode（弹性模式？）
当我们非线性对约束进行线性化时，可能会导致子问题没有可行域，即没有可行解(Infeasible)。例如，有两个约束分别为
$$
x\leq{1},\quad x^{2}\leq 4
$$
在$x_{k}=1$处线性化为：
$$
-p\geq 0 \quad 2p-3\geq 0 
$$
可知该子问题没有可行域，或无可行解。

我们可以采取类似$l_{1}$惩罚问题的方式重新构建该问题为
$$
\begin{align}
\min_{x,v,w,t} f(x)&+\mu \sum_{i\in\mathcal{E}}(v_{i}+w_{i})+\mu \sum_{i\in \mathcal{I}}t_{i} \\
s.t. \quad c_{i}(x) &= v_{i}-w_{i},\space i\in \mathcal{E}, \\
c_{i}(x)&\geq-t_{i}, \space i\in \mathcal{I}, \\
v,w,t &\geq 0,
\end{align}
$$
其中$\mu>0$为惩罚因子。该子问题始终是可行的(feasible)，软件SNOPT即采用此方法。

2）二阶矫正
有时候线性化的约束可能不够准确，将非线性约束展开到二阶进行矫正。并用$p^{T}_{k}\nabla^{2}c_{i}(x_{k})p_{k}$代替$p^{T}\nabla^{2}c_{i}(x_{k})p$ 。 

成熟软件：
- 线搜索法SQP：SNOPT
- 信赖域f法SQP：FILTERSQP
- SLQP：KNITRO

### 9 内点法(障碍法)

考虑问题
$$
\begin{align}
\min_{x, s} &f(x) \\
s.t.\quad c_{E}(x) &= 0, \\
c_{I}(x) - s &= 0, \\
s&\geq 0
\end{align}

$$
引入$\log$函数将原问题近似为

$$
\begin{align}
\min_{x, s} f(x)-&\mu \sum_{i=1}^{m}\log(s_{i}) \\
s.t.\quad c_{E}(x) &= 0, \\
c_{I}(x) - s &= 0,
\end{align}

$$
其中$\mu >0$.

该问题的KKT条件为

$$
\begin{align}
\nabla_{x}L = \nabla f(x) - A_{E}^{T}(x)y - A_{I}^{T}(x)z &= 0 \\
\nabla_{s}L = -\mu S^{-1}e + z &=0 \\
c_{E}(x) &= 0 \\
c_{I}(x)-s &=0
\end{align}


$$

迭代求解下降方向：
$$
\begin{align}
\begin{bmatrix}
\nabla_{xx}^{2}L & 0 & -A^{T}_{E}(x) & -A_{I}^{T}(x) \\
0 & Z & 0 & S \\
A_{E}(x) & 0 & 0 & 0 \\
A_{I}(x) & -I & 0 & 0 
\end{bmatrix}\begin{bmatrix}
p_{x} \\
p_{s} \\
p_{y} \\
p_{z}
\end{bmatrix}= \\

-\begin{bmatrix}
\nabla f(x)-A_{E}^{T}(x)y-A_{I}^{T}(x)z \\
Sz - \mu e \\
c_{E}(x) \\
c_{I}(x)-s
\end{bmatrix}

\end{align}
$$

这里的$L = f(x)-y^{T}c_{E}(x) - z^{T}(c_{I}(x)-s)$, $y, z$分别为等式约束和不等式约束的拉格朗日乘子。


### 10 MM 算法

![[Pasted image 20230718104039.png]]

![[Pasted image 20230717163033.png]]

算法思路：
原问题：
$$
\begin{align}
\min f(x) \\
s.t. \quad x\in \mathcal{H}
\end{align}
$$
其中$\mathcal{H}$是非空闭集，$f$是连续函，设当$\Vert{x}\Vert\to \infty$时，$f(x)\to \infty$。

迭代计算分为两步：

step1: 构造一个代理函数$g(x|x_{t})$, 满足以下条件：
$$
g(x|x_{t})\geq f(x)+c_{t},\forall x\in \mathcal{H}
$$
其中$c_{t} = g(x_{t}|x_{t})-f(x_{t})$。

step2: 
$$
x_{t+1} = \mathop{\arg \min_{x\in\mathcal{H}}\space g(x|x_{t})}
$$
序列$f(x_{t})_{t\in \mathbb{N}}$ 是非增的，因为：
$$
f(x_{t+1})\leq g(x_{t+1}|x_{t})-c_{t}\leq g(x_{t}|x_{t})-c_{t}=f(x_{t})
$$

应用MM思路的算法举例：SCA(Successive Convex Approximation) Algorithm.

1）近似性能指标
考虑问题
$$
\begin{align}
&\min_{x} \space f(x)+h(x) \\
&s.t. \space x\in \mathcal{X}
\end{align}
$$
其中$f:\mathcal{X}\to \mathbb{R}$的导数是Lipschitz连续的，$h:\mathcal{X}\to \mathbb{R}$是凸的。

代理函数$g(x|x_{t})$是凸的且满足$\nabla g(x|x_{t}) = \nabla f(x_{t})$, 近似的子问题为
$$
\begin{align}
&\min_{x}\space g(x|x_{t})+\frac{\tau}{2}(x-x_{t})^{T}Q(x_{t})(x-x_{t}
)+h(x) \\
&s.t. \space x\in \mathcal{X}
\end{align}

$$
其中$Q(x_{t})\in \mathbb{S}_{++}$

2）近似性能指标和约束
考虑问题
$$
\begin{align}
&\min_{x} \space f_{0}(x)
&s.t.\space f_{i}(x)\leq 0,i=1,\cdots,m
\end{align}
$$
设$f_{i}$是可导的，则在第$t$次迭代时，求解如下凸的子问题：
$$
\begin{align}
&\min_{x}\space g_{0}(x|x_{t}) \\
&s.t. \space g_{i}(x|x_{t})\leq0, i=1,\cdots,m, \\ 
\end{align}
$$
其中$g_{i}(x|x_{t}),i=0,\cdots,m$是凸函数且满足
$$
\begin{align}
g_{i}(x_{t}|x_{t}) &= f_{i}(x_{t}) \\
\nabla g_{i}(x_{t}|x_{t}) &= \nabla f_{i}(x_{t}) \\
g_{i}(x|x_{t})&\geq f_{i}(x) \\
\end{align}
$$
若序列$(x_{t})_{t\in \mathbb{N}}$收敛，则其极限为KKT点。



### 11 序列凸优化法

采用信赖域法框架

![[Pasted image 20230717155912.png]]


### 12 轨迹优化问题（大规模非线性规划问题）

1）EXAMPLE 1

$$
\begin{align}
\min J &=\int_{t_{0}}^{t_{f}} \frac{1}{2}u^{2} \, \mathrm{d}t  \\
s.t.\quad \dot{x} &= f(x,t) \\
c_{min}&\leq c \leq c_{max} \\
u_{min}&\leq u \leq u_{max} \\
x(t_{0})&=x_{0} \\
x(t_{f}) &= x_{f}
\end{align}
$$
设使用梯形积分进行转换，则NLP问题为：
$$
\begin{align}
\min J &= \sum_{k=0}^{N-1} \frac{h_{k}}{2}(u_{k}^{2}+u_{k+1}^{2}), \\
x&=[x_{0},x_{1},\cdots, x_{N}], \\
u&=[u_{0},u_{1},\cdots,u_{N}],  \\
s.t.\quad  x_{k+1}-x_{k}&= \frac{1}{2}h_{k}(f(x_{k},u_{k})+f(x_{k+1},u_{k+1}))\\
c_{min}&\leq c(x,u) \leq c_{max} \\
u_{min}&\leq u \leq u_{max} \\

\end{align}
$$

其中$N$为划分的子区间数，即将时间划分为$t_{0},t_{1},\cdots,t_{N}$ , $h_{k} = t_{k+1}-t_{k}$。
