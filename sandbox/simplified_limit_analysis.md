# Simplified Limit Analysis for Dirichlet-Multinomial Log-PMF

## Reparametrized Log-PMF with Clear Limit Structure

Starting from the exact log-PMF and applying the Beta function approximation:

\begin{equation}\label{eq:reparametrized}
\log P(\mathbf{y} | \tau, N) = \log\left(\frac{N!}{K^N \prod_{k=1}^K y_k!}\right) + \underbrace{\log \Gamma(1/\tau) - \log \Gamma(N + 1/\tau)}_{A(\tau)} + \underbrace{\sum_{k=1}^K \left[\log \Gamma(y_k + 1/(K\tau)) - \log \Gamma(1/(K\tau))\right]}_{B_k(\tau)}
\end{equation}

## Asymptotic Expansion for Small $\tau$

### Term A: $\log \Gamma(1/\tau) - \log \Gamma(N + 1/\tau)$

Using Stirling's approximation for large $z = 1/\tau$:

\begin{equation}\label{eq:stirling_diff}
\log \Gamma(z) - \log \Gamma(z + N) \approx -N \log z + O(1/z)
\end{equation}

Thus:

\begin{equation}\label{eq:termA}
A(\tau) = \log \Gamma(1/\tau) - \log \Gamma(N + 1/\tau) \approx -N \log(1/\tau) + O(\tau) = N \log \tau + O(\tau)
\end{equation}

### Term B: $\log \Gamma(y_k + 1/(K\tau)) - \log \Gamma(1/(K\tau))$

Similarly:

\begin{equation}\label{eq:log_gamma_diff}
\log \Gamma(z + y_k) - \log \Gamma(z) \approx y_k \log z + O(1/z) \quad \text{where} \quad z = 1/(K\tau)
\end{equation}

Thus:

\begin{equation}\label{eq:termB}
B_k(\tau) = \log \Gamma(y_k + 1/(K\tau)) - \log \Gamma(1/(K\tau)) \approx y_k \log(1/(K\tau)) + O(\tau) = y_k \log(1/\tau) - y_k \log K + O(\tau)
\end{equation}

### Combined Asymptotic Behavior

Summing all terms:

\begin{multline}\label{eq:combined}
\log P(\mathbf{y} | \tau, N) \approx \log\left(\frac{N!}{K^N \prod_{k=1}^K y_k!}\right) + N \log \tau \\
+ \sum_{k=1}^K \left(y_k \log(1/\tau) - y_k \log K\right) + O(\tau)
\end{multline}

Simplifying the $\log \tau$ terms:

\begin{multline}\label{eq:simplified}
= \log\left(\frac{N!}{K^N \prod_{k=1}^K y_k!}\right) + N \log \tau + \left(\sum_{k=1}^K y_k\right) \log(1/\tau) - \left(\sum_{k=1}^K y_k\right) \log K + O(\tau)
\end{multline}

Since $\sum_{k=1}^K y_k = N$:

\begin{equation}\label{eq:cancellation}
= \log\left(\frac{N!}{K^N \prod_{k=1}^K y_k!}\right) + N \log \tau + N \log(1/\tau) - N \log K + O(\tau)
\end{equation}

The $\log \tau$ terms cancel exactly:

\begin{equation}\label{eq:after_cancel}
= \log\left(\frac{N!}{K^N \prod_{k=1}^K y_k!}\right) - N \log K + O(\tau)
\end{equation}

## Final Simplified Form

\begin{equation}\label{eq:final}
\log P(\mathbf{y} | \tau, N) = \log\left(\frac{N!}{K^N \prod_{k=1}^K y_k!}\right) - N \log K + O(\tau)
\end{equation}

## Limit as $\tau \to 0$

### Leading Order Behavior

As $\tau \to 0^+$:

\begin{equation}\label{eq:limit}
\log P(\mathbf{y} | \tau, N) \to \log\left(\frac{N!}{K^N \prod_{k=1}^K y_k!}\right) - N \log K
\end{equation}

### Exact Multinomial PMF

Exponentiating gives the uniform multinomial distribution:

\begin{equation}\label{eq:multinomial}
P(\mathbf{y} | \tau = 0, N) = \frac{N!}{K^N \prod_{k=1}^K y_k!}
\end{equation}

## Interpretation

This simplified analysis clearly shows:

1. **Divergent Terms**: The $N \log \tau$ and $N \log(1/\tau)$ terms cancel exactly
2. **Finite Limit**: Only the multinomial term remains as $\tau \to 0$
3. **Correction Terms**: The $O(\tau)$ terms vanish in the limit
4. **Convergence**: The Dirichlet-multinomial $\to$ multinomial as overdispersion vanishes

The key insight is that the apparently divergent logarithmic terms in $\tau$ cancel precisely, leaving a well-defined multinomial distribution in the $\tau \to 0$ limit.