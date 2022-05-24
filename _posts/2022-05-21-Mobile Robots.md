---
layout:     post
title:      "移动机器人定位问题"
# subtitle:   " "\概率，概率，还是概率\""
date:       2022-05-21 10:45:00
author:     "Lu_fm"
header-img: "img/sr_dog.png"
catalog: true
tags:
    - Control Science
    - Robotics
---


### 定位问题描述
$
p(x_t|Z^t,U^{t-1},m) = \eta p(z_t|x_t,m)\int p(x_t|x_{t-1},u_{t-1},m)p(x_{t-1}|Z^{t-1}U^{t-2},m)dx_{t-1}
$

---



### 运动模型

$$
p(x_t|x_{t-1},u_{t-1},m) = \eta p(m|x_t,x_{t-1},u_{t-1})p(x_t|x_{t-1},u_{t-1})
$$

- 右式第一项为机器人在控制指令$u_{t-1}$从$x_{t-1}$移动到$x_t$的条件下获得地图的可能性，如果发生碰撞则可能性低，此项可降低运动模型的不确定性,如果不考虑此项，不确定性会无上限不断增大。

  - 计算复杂，无法实际应用时可简化为：

  	$$
    p(x_t|x_{t-1},u_{t-1},m) \approx(x_t|x_{t-1},u_{t-1})\\
    
    
    p(x_t|x_{t-1},u_{t-1},m) = \eta p(x_t|m)p(x_t|x_{t-1},u_{t-1})
    $$

    

- 第二项为根据上一时刻的指令和位置估计当前位姿的概率分布。若为速度运动模型，则$u_{t-1}$为发送给机器人的指令，若为里程计模型，则$u_{t-1}$为里程计数据。

  - 里程计运动模型（航位推算）

	$$
    \begin{bmatrix}
    x'\\
    y'\\
    \theta'
    \end{bmatrix}=\begin{bmatrix}
    x+\delta_{trans}cos(\theta+\delta_{rot1})\\
    y+\delta_{trans}sin(\theta+\delta_{rot1})\\
    \theta+\delta_{rot1}+\delta_{rot2}
    \end{bmatrix}
  $$
    真实值 = 测量值 - 误差， 以$\delta_{rot1}$为例，其余类似
    $
    \hat \delta_{rot1} = \delta_{rot1}-\tilde\delta_{rot1}\\
    $
    闭式求解：
  
    1. $u_{t-1}$ 表示为$\delta=(\delta_{rot1},\delta_{trans},\delta_{rot2})^T$
    
    2. 根据$x_t, x_{t-1}$计算实际的运动变化量$\hat\delta=(\hat\delta_{rot1},\hat\delta_{trans},\hat\delta_{rot2})^T$
    
    3. 根据误差分布估计$x_t$ 的分布：
       $
       p(x_t|,x_{t-1},u_{t-1}) =  p(\tilde\delta_{rot1})p(\tilde\delta_{trans})p(\tilde\delta_{rot2})
       $
       假设噪声为高斯分布，即有：
       
  $$
       \begin{align}
       &p(\tilde\delta_{rot1}) = N(\tilde\delta_{rot1};0,\alpha_1|\delta_{rot1}|+\alpha_2|\delta_{trans}|)\\
       &p(\tilde\delta_{trans}) = N(\tilde\delta_{trans};0,\alpha_3|\delta_{trans}|+\alpha_4|\delta_{rot1}+\delta_{rot2}|)\\
       &p(\tilde\delta_{rot2}) = N(\tilde\delta_{rot2};0,\alpha_5|\delta_{rot2}|+\alpha_6|\delta_{trans}|)
       \end{align}
  $$
  
   共有6个参数需要选择。
  
   4.随机采样求解，将噪声建模为高斯分布。
  
   - 在噪声的高斯分布中随机采样, 获得控制指令

  $$
  \tilde\delta_{rot1}=sample(\alpha_1|\delta_{rot1}|+\alpha_2|\delta_{trans}|)\cdots
  $$
  
   - 利用上述里程运动模型计算$x_t$ 分布

   - 重复上述过程，所得$x_t$点集合构成对$p(x_t \|x_{t-1},u_{t-1})$的描述。
  


### 观测模型$p(z_t|x_t,m)$

观测模型依赖于传感器模型。

- 特征的传感器模型
  - 特征地图表示为$m=\{m_1,m_2,\cdots,m_N\}$

​			观测$z_t = \{f_t^1,f_t^2,\cdots,f_t^K\}$, $f_t^i \triangleq (r_t^i,\varphi_t^i,s_t^i)$ , $s_t^i$为特征标识。
$
p(z_t|x_t,m) =p(f_t^1,f_t^2,\cdots,f_t^K|x_t,m) = \Pi_{i=1}^Kp(f_t^i|x_t,m)
$
​将地图特征$m_j$转换到机器人坐标系下即为观测特征真值。

$$
    \hat f_t^i = 
    \begin{bmatrix}
    \sqrt{(m_{jx}-x)^2)+(m_{jy}-y)^2}\\
    atan2(m_{jy}-y, m_{jx}-x)\\
    s_j
    \end{bmatrix}\\
    
    f_t^i = \hat f_t^i + \tilde \delta_f
$$


根据观测误差模型，则

  $$ p(f_t^i|x_t,m) = N(f_t^i;\hat f_t^i,\sigma_f^2) $$

- 基于物理建模的激光传感器模型

  - $$
    p(z_t^k|x_t,m) = \alpha_{hit}p_{hit}(z_t^k|x_t,m)+\alpha_{short}p_{short}(z_t^k|x_t,m)+\alpha_{max}p_{max}(z_t^k|x_t,m)+\alpha_{rand}p_{rand}(z_t^k|x_t,m)
    $$

    

归一化
	$
\alpha_{hit}+\alpha_{short}+\alpha_{max}+\alpha_{rand}=1
	$
	
hit： 带少量噪声的正确测量\\
short: 临时障碍的测量\\
max: 达到最大距离的测量\\
rand: 随机错误测量



### EKF自定位法

#### 1.KF

- 运动方程

	$$
x_t =A_tx_{t-1}+B_tu_{t-1}+v_t\\
x_{t-1}\in R^n,  x_{t-1}\sim N(\mu_{t-1},\Sigma_{t-1})\\
u_{t-1}\in R^m\\
v_t\in R^n, v_t \sim N(0,R_t)\\
A_t\in R^{n\times n}, B_t\in R^{n\times m}\\
p(x_t|x_{t-1},u_{t-1}) = N(x_t;A_tx_{t-1}+B_t u_{t-1},R_t)
	$$

- 观测方程

  $$
  z_t = C_tx_t + \omega_t,\quad C_t\in R^{k\times n}, \omega_t\sim N(0,Q_t)\\
  p(z_t|x_t) = N(z_t|C_tx_t,Q_t)
  $$

- 状态预估

  $$
  \overline p(x_t|Z^{t-1},U^{t-1}) = N(x_t;\overline \mu_t,\overline \Sigma_t)\\
  \overline \mu_t = A_t\mu_{t-1}+B_t u_{t-1}\\
  \overline \Sigma_t = A_t\Sigma_{t-1}A_t^T+R_t
  $$
  
- 观测更新

  $$
  p(x_t|Z^t,U^{t-1}) = \eta p(z_t|x_t)\overline p(x_t|Z^{t-1},U^{t-1}) \\
  
   p(x_t|Z^t,U^{t-1}) = N(x_t|\mu_t,\Sigma_t)\\
  
   \mu_t = \overline \mu_t+K_t(z_t-C_t\overline \mu_t)\\
   \Sigma_t = [I-K_tC_t]\overline \Sigma_t\\
  
  K_t = \overline \Sigma_tC_t^T(C_t\overline \Sigma_t^T+Q_t)^{-1},新信息矩阵
  $$

局限：要求线性函数

---



#### 2.EKF

采用非线性函数表示运动方程和观测方程，并计算非线性转换后实际概率的高斯近似。

$$
    x_t = g(x_{t-1},u_{t-1})+v_t\\
    z_t = h(x_t)+\omega_t
$$


用$x_{t-1}$处的一阶泰勒展开近似：

$$
    x_t =\approx g(\mu_{t-1},u_{t-1})+(x_{t-1},u_{t-1})G_t+v_t\\
    G_t = \frac{\partial g(x_{t-1},u_{t-1})}{\partial x_{t-1}}|x_{t-     1}=\mu_{t-1}
$$

观测方程也为：

$$
    h(x_t)\approx h_(\overline \mu_t)+(x_t-\overline \mu_t)H_t\\
    H_t = \frac{\partial h(x_t)}{\partial x_t}|x_t = \overline \mu_t =   g(\mu_{t-1},u_{t-1})
$$

位姿预估

$$
    \overline p(x_t|Z^{t-1},U^{t-1}) = N(x_t;\overline \mu_t,\overline 		\Sigma_t)\\
	\overline \mu_t = g(\mu_{t-1},u_{t-1})\\
	\overline \Sigma_t = G_t\Sigma_{t-1}G_t^T+R_t
$$

观测更新：

$$
	p(x_t|Z^t,U^{t-1}) = \eta p(z_t|x_t)\overline p(x_t|Z^{t-1},U^{t-1})\\

 	p(x_t|Z^t,U^{t-1}) = N(x_t|\mu_t,\Sigma_t)\\

	\mu_t = \overline \mu_t+K_t(z_t-h(\overline \mu_t))\\

	\Sigma_t = [I-K_tH_t]\overline \Sigma_t\\

	K_t = \overline \Sigma_t H_t^T(H_t\overline\Sigma_tH_t^T+Q_t)^{-1}, 新信息矩阵
$$
