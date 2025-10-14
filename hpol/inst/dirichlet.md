---
title: cc multinom dirichelet
---



Strata
$Y_{i\ell} = \{Y_{i\ell m}, m = 1 \ldots M_{i\ell} \}$ 

$$
\begin{aligned}
\lambda_{ij}(t) = & Z_i(t) \rho_{ij}(t) \mu_i(t)  \\
Z_{i}(t) \sim &  \text{Gamma}(1/\nu^2, 1/\nu^2)\\
Y_{i\ell m} \sim & \text{Poisson}[\lambda_{ij}(t_{i\ell m})]\\
Y_{i\ell} | \text{sum}\{Y_{i\ell}\},  Z, \mu \sim & 
\text{Multinom}\left(\text{sum}\{Y_{i\ell}\}, P_{i\ell}\right)\\
P_{i\ell m} | Z, \eta \propto & Z_{i}(t_{i\ell m})\mu_i(t_{i\ell m})\\
\end{aligned}
$$
Note $\text{E}[Z_i(t)] = 1$, $\text{sd}[Z_i(t)] = \nu$ and
$$
\begin{aligned} 
P_{i\ell} | \eta \sim & \text{Dirichelet}\left[1/\nu^2, \bar\mu_{i\ell}\right] \\
\bar \mu_{i\ell m}  = & \mu_{i}(t_{i\ell m}) \left/ \sum_m \mu_{i}(t_{i\ell m}) \right. \\
\end{aligned}
$$


Dirichlet multinomial 
(\href{https://w.wiki/D$by}{wikipedia})
$$
\begin{aligned} 
Y_{i\ell} | \text{sum}\{Y_{i\ell}\} , \mu \sim & \text{DirMult}\left[\text{sum}\{Y_{i\ell}\}, 1/\nu^2 , \bar\mu_{i\ell} \right]\\
\pi(Y_{i\ell}; 1/\nu^2, \bar\mu_{i\ell}) = & 
\frac{
    \Gamma(1+\text{sum}\{Y_{i\ell}\}  )\Gamma(1/\nu)
}{
    \Gamma(1/\nu^2+\text{sum}\{Y_{i\ell}\}  )
} 
    \prod_{m=1}^{M_{i\ell}} \frac{\Gamma(Y_{i\ell m}+\bar\mu_{i\ell m}/\nu^2)
}{
    \Gamma(Y_{i\ell m}+1)\Gamma(\bar\mu_{i\ell m}/\nu^2)
} \\
= & \text{sum} \{Y_{i\ell}\}  B\left(1/\nu^2, \text{sum} \{Y_{i\ell}\} \right)
\left/
\prod_{m; Y_{i\ell m}>0} Y_{i\ell m}B\left(\bar\mu_{i\ell m}/\nu^2, Y_{i\ell m}\right)
\right.
\end{aligned}
$$

Log likelikhood
$$
\begin{aligned}
\log \pi(Y_{i\ell}; 1/\nu^2, \bar\eta_{i\ell}) = &
\log \Gamma \left(1+\text{sum} \{Y_{i\ell}\}\right) + \log \Gamma(1/\nu^2) 
- \log  \Gamma(1/\nu^2+\text{sum}\{Y_{i\ell}\}  ) \\
& + \sum_m  \log \Gamma(Y_{i\ell m}+\bar\mu_{i\ell m}/\nu^2) 
- \log \Gamma(Y_{i\ell m}+1) -\log \Gamma(\bar\mu_{i\ell m}/\nu^2)
\end{aligned}
$$
