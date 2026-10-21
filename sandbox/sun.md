
$$
\begin{aligned}
\begin{bmatrix}
\Sigma & \Delta^T \\
\Delta & \Omega
\end{bmatrix}
= &
\begin{bmatrix}
I & 0 \\
L & I
\end{bmatrix}
\begin{bmatrix}
A_1 & 0 \\
0 & A_2
\end{bmatrix}
\begin{bmatrix}
I & L^T \\
0 & I
\end{bmatrix}\\
= & 
\begin{bmatrix}
A_1 & A_1 L_1^T \\
L_1 A_1 &  L_1 A_1 L_1^T + A_2
\end{bmatrix}
\end{aligned}
$$

- Want $var(U|V>0) = I$, mean = 0
- what conditions on A1, A2 and L?
- choose some valid A1, A2, L, and simultate and scatterplot
- find some cases where mode is a long way from mean