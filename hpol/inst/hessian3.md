---
title: deriv of laplace approx
---


[TMB](https://www.jstatsoft.org/article/view/v070i05) page 8

[matrix cookbook](https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf) section 2

Derivatives and tensor
$$
\begin{aligned}
T_{\theta UU} =     & \frac{\partial }{\partial U} H(U, \theta)\bigg|_{U = \hat U(\theta)} =
  \frac{\partial^3}{(\partial U )^2\partial\theta} \log \pi[ Y| U;\theta] \bigg|_{U = \hat U(\theta)} \\
H_{UU} =     & H\left(\hat U(\theta), \theta\right) = 
 \frac{\partial^2}{(\partial U )^2} \log \pi[ Y| U;\theta] \bigg|_{U = \hat U(\theta)}\\
H_{U\theta} =     &  \frac{\partial^2}{\partial U \partial \theta} \log \pi[ Y| U;\theta] \bigg|_{U = \hat U(\theta)}
\end{aligned}
$$


$$
\begin{aligned}
\partial \log |X| = & \text{Tr}(X^{-1} \partial X) \\
\partial \log |H_{UU}[\hat U(\theta), \theta]  / \partial \theta = & 
\end{aligned}
$$

Random effects

TMB says
$$
\begin{aligned}
\frac{\partial}{\partial  \theta} \hat U(\theta) = & - H\left(\hat U(\theta), \theta\right)^{-1}    \cdot 
    \frac{\partial^2}{\partial U \partial \theta} \log \pi[ Y| U;\theta] \bigg|_{U = \hat U(\theta)}\\
    = & - H_{UU}^{-1} H_{U \theta}
\end{aligned}
$$

Consider

$$
\begin{aligned}
 U^{(1)}(\theta, U^{(0)}) = & U^{(0)} + H\left( U^{(0)}, \theta\right)^{-1}  G\left( U^{(0)}, \theta\right) \\
 \frac{\partial}{\partial  \theta} \hat U^{(1)}(\theta, U^{(0)}) = & 
 \frac{\partial}{\partial  \theta} \left\{H\left( U^{(0)}, \theta\right)^{-1} \right\}
    G\left( U^{(0)}, \theta\right) + 
  H\left( U^{(0)}, \theta\right)^{-1} H^{(0)}_{U\theta} \\
 \frac{\partial}{\partial  \theta_p} \hat U^{(1)}(\theta, U^{(0)})   = & - H^{-1}  T_{\cdot \cdot p} H^{-1} 
    G\left( U^{(0)}, \theta_p\right) + 
  H^{-1} H_{U\theta_p} \\
  = & - H^{-1} \left[ T_{\cdot \cdot p}  H^{-1} 
    G\left( U^{(0)}, \theta_p\right)  + H_{U\theta_p} \right]
\end{aligned}
$$

$G=0$ and we can ignore


Hessian


$$
\begin{aligned}
  \frac{\partial}{\partial \theta_p}  H_{ij}\left(\hat U(\theta), \theta\right) \   = & 
     \frac{\partial}{\partial \theta_p} H_{ij}(U, \theta)|_{U = \hat U(\theta)} +
      \frac{\partial}{\partial U} H_{ij}(U, \theta)|_{U = \hat U(\theta)} \frac{\partial}{\partial \theta_p}
      \hat U(\theta) \\
      = &  T_{ijp} + \sum_k T_{ijk} 
      \hat U'_{kp}
\end{aligned}
$$

compute $\hat U'$, 

log det:


$$
\begin{aligned}
  \frac{\partial}{\partial \theta_p} \log \left| H\left(\hat U(\theta), \theta\right) \right|   = & 
    \text{trace}\left\{
     H\left(\hat U(\theta), \theta\right)^{-1}  
     \frac{\partial }{\partial \theta_p} H\left[ \hat U(\theta), \theta\right] 
     \right\}\\
= & \sum_{ij} (H^{-1})_{ij}   \frac{\partial}{\partial \theta_p}  H_{ij}\left(\hat U(\theta), \theta\right) \\
= & \sum_{ij} (H^{-1})_{ij} T_{ijp} + \sum_{ijk} (H^{-1})_{ij} T_{ijk}    \hat U'_{kp}
\end{aligned}
$$


TMB's log det

$$
\begin{aligned}
  \frac{\partial}{\partial \theta_p} \log \left| H\left(\hat U(\theta), \theta\right) \right|   = & 
    \text{trace}\left\{
     H\left(\hat U(\theta), \theta\right)^{-1}  
     \frac{\partial }{\partial \theta_p} H\left( U, \theta\right) \bigg|_{U = \hat U(\theta)}
     \right\}\\
     = &    \text{trace}\left\{
     H_{UU}^{-1} \cdot  T_{\theta_p U U} \right\} 
\end{aligned}
$$

That ignores changes in $\hat U(\theta)$.  Consider

$$
\begin{aligned}
  \frac{\partial}{\partial \theta_p} \log \left| H\left(\hat U(\theta), \theta\right) \right|  = &
    \frac{\partial}{\partial U} 
        \log \left| H\left(U, \theta\right) \right| \bigg|_{U = \hat U(\theta)} 
            \frac{\partial}{\partial  \theta_p} \hat U(\theta)  + 
  \frac{\partial}{\partial \theta_p} 
    \log \left| H\left( U, \theta\right) \right| \bigg|_{U = \hat U(\theta)} \\
\frac{\partial}{\partial U_i} 
        \log \left| H\left(U, \theta\right) \right| \bigg|_{U = \hat U(\theta)} = & \text{trace}\{
       H_{UU}^{-1} \cdot T_{UUU_i} \} \\
         \frac{\partial}{\partial \theta_p} 
    \log \left| H\left(U, \theta\right) \right| \bigg|_{U = \hat U(\theta)} =&
     \text{trace}\{ H_{UU}^{-1} \cdot T_{U U \theta_p }\}\\
       \frac{\partial}{\partial \theta_p} \log \left| H\left(\hat U(\theta), \theta\right) \right|  = & \sum_i  \text{trace}\{
       H_{UU}^{-1} \cdot T_{UUU_i} \}  e_i H_{UU}^{-1} H_{U\theta_p} +      \text{trace}\{ H_{UU}^{-1} \cdot T_{U U \theta_p }\} \\
    %   = & V^T H_{UU}^{-1} H_{U\theta_p} + \text{trace}\{ H_{UU}^{-1} \cdot T_{U U \theta_p }\}\\
       V_i = &   \text{trace}\{
       H_{UU}^{-1} \cdot T_{UUU_i} \} \\
       W_p = &  \text{trace}\{ H_{UU}^{-1} \cdot T_{U U \theta_p }\}\\
\frac{\partial}{\partial \theta} \log \left| H\left(\hat U(\theta), \theta\right) \right| = &
V H_{UU}^{-1} H_{U\theta} + W
\end{aligned}
$$

Data component
$$
\begin{aligned}
 \frac{\partial}{\partial \theta} \log \pi[ Y|\hat U(\theta);\theta] = &
 \frac{\partial}{\partial U} \log \pi[ Y| U;\theta] \bigg|_{U = \hat U(\theta)} 
 \frac{\partial}{\partial  \theta} \hat U(\theta) + 
 \frac{\partial}{\partial \theta} \log \pi[ Y| U ;\theta] \bigg|_{U = \hat U(\theta)} \\
= & - G_U H_{UU}^{-1} H_{U\theta} + G_\theta
\end{aligned}
$$


Likelihood
$$
\begin{aligned}
\ell(\theta) = &\log \pi[ Y|\hat U(\theta);\theta] + 0.5 \log | H(\hat U(\theta), \theta) | \\
\frac{\partial}{\partial \theta} \ell(\theta) = & \frac{\partial}{\partial \theta} \log \pi[ Y|\hat U(\theta);\theta] + 
    0.5 \frac{\partial}{\partial  \theta}  \log | H\left(\hat U(\theta), \theta\right) | \\
    = & - G_U H_{UU}^{-1} H_{U\theta} + G_\theta + V H_{UU}^{-1} H_{U\theta} + W \\
    = & G_\theta + V H_{UU}^{-1} H_{U\theta} + W
\end{aligned}
$$

TMB ignores $V H_{UU}^{-1} H_{U\theta}$



# fourth

$$
\begin{aligned}
\frac{\partial^2}{\partial \theta_q\,\partial \theta_p}
\log\big|H(\hat U(\theta),\theta)\big|
&=
\sum_{a,b}
\frac{\partial \hat U_a}{\partial \theta_q}
\frac{\partial \hat U_b}{\partial \theta_p}
\operatorname{tr}\!\big(H^{-1} F_{ab..}\big)
   \operatorname{tr}\!\big(H^{-1} T_{b..} H^{-1} T_{a..}\big)
\\
& + \quad
   \sum_a
\frac{\partial \hat U_a}{\partial \theta_p}
\operatorname{tr}\!\big(H^{-1} F_{aq..}\big)
  \operatorname{tr}\!\big(H^{-1} T_{q..} H^{-1} T_{a..}\big)
\\ 
& + \quad
   \sum_a
\frac{\partial \hat U_a}{\partial \theta_q}
\operatorname{tr}\!\big(H^{-1} F_{ap..}\big)
   \operatorname{tr}\!\big(H^{-1} T_{p..} H^{-1} T_{a..}\big)
\\ 
& + \quad
   \sum_a
\frac{\partial^2 \hat U_a}{\partial \theta_q\,\partial \theta_p}
\,\operatorname{tr}\!\big(H^{-1} T_{a..}\big)
\\ 
& + \quad
   \operatorname{tr}\!\big(H^{-1} F_{pq..}\big)
   \operatorname{tr}\!\big(H^{-1} T_{q..} H^{-1} T_{p..}\big)
\end{aligned}
$$


$$
\begin{aligned}
\frac{\partial^2 \hat U_a}{\partial \theta_q\,\partial \theta_p}
&= - \sum_b (H^{-1})_{ab} \Bigg[
\sum_{c,d} T_{bcd} \frac{\partial \hat U_c}{\partial \theta_q}
\frac{\partial \hat U_d}{\partial \theta_p} + 
\sum_{c} T_{bcq} \frac{\partial \hat U_c}{\partial \theta_p} + 
\sum_{c} T_{bcp}
\frac{\partial \hat U_c}{\partial \theta_q} +  T_{bpq} \Bigg]
\end{aligned}
$$