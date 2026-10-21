# Concise Derivation of Reparametrized Dirichlet-Multinomial Log-PMF

## Beta Function Infinite Product Representation

The Beta function has an infinite product representation:

\begin{equation}
B(x, y) = \frac{x + y}{xy} \prod_{n=1}^\infty \left(1 + \frac{xy}{n(x + y + n)}\right)^{-1}
\end{equation}

### Logarithm of Beta Function

Taking logarithms gives:

\begin{equation}
\log B(x, y) = \log(x + y) - \log x - \log y - \sum_{n=1}^\infty \log\left(1 + \frac{xy}{n(x + y + n)}\right)
\end{equation}

### Special Case: $\log B(a, 1/b)$

Substituting $y = 1/b$ and simplifying:

\begin{equation}
\log B\left(a, \frac{1}{b}\right) = \log(ab + 1) - \log a - \sum_{n=1}^\infty \log\left(1 + \frac{a}{n(ab + 1 + bn)}\right)
\end{equation}

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
