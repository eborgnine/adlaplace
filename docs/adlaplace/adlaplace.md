---
title: deriv of laplace approx
---

# Intro


Ad gives derivatives of $\log\pi(Y, U;\theta)$.  We need $\log\pi(Y;\theta)$.

Laplace approximation $h(\theta) = \log |H(\hat U(\theta), \theta)|$ deriv of $h(\theta)$ is a third order derivative.  

Goal: do it in a way that's compatible with special functions.  clear separation of general-purpose backend and minimal frontend, suitable for modification and extension.  modular, combine with aghq, rstan.

## Motivating example

Neg binom glmm, log gamma function.  know derivatives, psi gamma


# Preliminaries

## AD

- chain rule, tape
- Atomics, code analytitcal derivatives i.e. cosines, reduce tape size
- Reverse steps, contractions: compute $A^T H$, specify $A = e_p$ for $p$ th row of H
- Sparsity, algorithms to get $H$ for a specified sparsity pattern

## GLMM

$\hat U(\theta)$ no analytical expression.  AD of N-R optimization

Identities. $\partial \log |X| = \text{Tr}(X^{-1} \partial X)|$

Double taping.  Atomics need to be compatible

# Methods

1. Decompose $X^{-1} = B B^T$. $\partial \log |X| = \text{Tr}(B \partial X B^T)|$, efficient using contractions
2. Parallelization, sparse shards of data, 


## Inner optimization $\hat U(\theta)$

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
  = & - H_{UU}^{-1} H_{Up}
\end{aligned}
$$

$G=0$ and we can ignore


# Determinant $\log|H(\hat U(\theta), \theta)|$


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
       H_{UU}^{-1} \cdot T_{UUi} \} \\
         \frac{\partial}{\partial \theta_p} 
    \log \left| H\left(U, \theta\right) \right| \bigg|_{U = \hat U(\theta)} =&
     \text{trace}\{ H_{UU}^{-1} \cdot T_{U U p }\}\\
       \frac{\partial}{\partial \theta_p} \log \left| H\left(\hat U(\theta), \theta\right) \right|  = & \sum_i  \text{trace}\{
       H_{UU}^{-1} \cdot T_{UUi} \}  e_i H_{UU}^{-1} H_{Up} +      \text{trace}\{ H_{UU}^{-1} \cdot T_{U U p }\} \\
       V_i = &   \text{trace}\{
       H_{UU}^{-1} \cdot T_{UUU_i} \} \\
       W_p = &  \text{trace}\{ H_{UU}^{-1} \cdot T_{U U \theta_p }\}\\
\frac{\partial}{\partial \theta} \log \left| H\left(\hat U(\theta), \theta\right) \right| = &
V H_{UU}^{-1} H_{U\theta} + W
\end{aligned}
$$

# Contractions

$H_{UU}^{-1} = B B^T$ 

$$
\begin{aligned}
\text{trace}\{H_{UU}^{-1} T_{UU\cdot}\} = & \text{trace}\{B T_{UU\cdot}B^T\} \\
 = & \sum_i B_i^T T_{UU\cdot} B_i\\
\end{aligned}
$$


# Parallelization and sparsity

$$
\log \pi(Y, U; \theta) = \sum_\ell \log \pi (Y^{(\ell)}| U; \theta) + \sum_m \log \pi(U^{(m)};\theta)
$$



$$
\text{trace}\{H_{UU}^{-1} T_{UU\cdot}\}  =  \sum_{i\ell} B^{ T}_i T^{(\ell)}_{UU\cdot} B_i + \sum_{im} B^{ T}_i T^{(m)}_{UU\cdot} B_i\\
$$

# Discussion

- Deriv of likelihood using third derivatives, compatible with atomics
- Efficient use of contractions, avoid direct computation of third tensor
- combination of sparsity and parallelization 
- common framework, clear back end




# References


[TMB](https://www.jstatsoft.org/article/view/v070i05) page 8

[matrix cookbook](https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf) section 2



# Appendix



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


# Code

AD objects

- GroupPack: for one shard ad function, work and sparsity, scatter maps.  cppp, no Rcpp.
- BackendContext: vector of GroupPack, integer vectors of hessian pattern, sizes.
- adpack api: simple c structure with exposed functions to call ad and return hessian pattern
- adpack handle: pointer to adpack api, accessible across packages

Creating the GroupPack's

- package's local function objectiveFunction.cpp that evaluates density functions
  - optionally sourcing logDensRandom.hpp and an empty logDensExtraEmpty.hpp from adLaplace
- sources adfun_create.hpp which 

