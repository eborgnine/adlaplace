### Where `U` and `V` come from

```cpp
f.Forward(0 , x0);      // 0-th-order coefficients  u(0)  ←  x0
f.Forward(1 , U);       // 1-st-order coefficients  u(1)  ←  U
f.Forward(2 , V);       // 2-nd-order coefficients  u(2)  ←  V
```

Internally CppAD now thinks you are studying the curve

$$
x(t)=x_0 \;+\; U\,t \;+\; V\,t^{2}.
$$

### What `Reverse(3 , w)` returns

For every independent variable $x_j$ it gives three numbers

$$
\begin{aligned}
dw[j\!*\!3+0]&=\frac{\partial W}{\partial u_j^{(0)}}\;,\\
dw[j\!*\!3+1]&=\frac{\partial W}{\partial u_j^{(1)}}\;,\\
dw[j\!*\!3+2]&=\frac{\partial W}{\partial u_j^{(2)}}\;,
\end{aligned}
$$

where

$$
W(u)=\sum_{k=0}^{2} w^{(k)}\;
     \frac{1}{k!}\,\frac{\partial^{\,k}}{\partial t^{k}}
     F\!\bigl(x(t)\bigr)\Big|_{t=0},
\qquad
u^{(k)}=\bigl\{u_j^{(k)}\bigr\}.
$$

---

## How `U` and `V` enter those formulas (scalar output)

Write $f_i   =\partial_{x_i}f,\;
      f_{ij} =\partial_{x_i x_j}^2f,\;
      f_{ijk}=\partial_{x_i x_j x_k}^3f$ at $x_0$.
Taylor-expand $f(x(t))$:

$$
\begin{aligned}
f(x(t)) &= f_0
        + \sum_{i}f_i   U_i\,t
        + \frac12\Bigl[\sum_{i}f_i V_i+\sum_{ij}f_{ij}U_iU_j\Bigr] t^{2}
        + O(t^{3}).
\end{aligned}
$$

Hence

$$
\begin{aligned}
y^{(0)} &= f_0,\\
y^{(1)} &= \sum_i f_i U_i,\\
y^{(2)} &= \sum_i f_i V_i+\sum_{ij}f_{ij}U_iU_j.
\end{aligned}
$$

Plugging these into $W(u)$ and differentiating gives

| returned component         | closed form (scalar case)                                                                                           |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **order 0**<br>`dw[j*3+0]` | $w^{(0)}f_j + w^{(1)}\sum_k f_{jk}\,U_k +\dfrac{w^{(2)}}{2}\Bigl(f_{j k}V_k+ \sum_{k\ell}f_{jk\ell}U_kU_\ell\Bigr)$ |
| **order 1**<br>`dw[j*3+1]` | $w^{(1)}f_j + w^{(2)}\sum_k f_{jk}\,U_k$                                                                            |
| **order 2**<br>`dw[j*3+2]` | $\dfrac{w^{(2)}}{2}\,f_j$                                                                                           |

---

### Key observations

* **`U` appears wherever a first-order coefficient is needed.**
  • If `w` selects $y^{(1)}$ or $y^{(2)}$, the $\,f_{jk}\,U_k$ and $f_{jk\ell}U_kU_\ell$ terms show up.
  • If `U` has a non-zero component only in $x_0$, every variable $x_j$ still “sees” that $U_0$ through the mixed derivatives $f_{j0}$ or $f_{jk0}$.

* **`V` can influence only the `dw[ * 0 ]` row** (because $y^{(2)}$ is the only place $V$ enters).

* **If your analytic Hessian/third-tensor is diagonal** but you still see cross-variable numbers, then at least one of:

  1. The tape really is not diagonal (possibly hidden couplings in your code).
  2. Some unintended non-zero component crept into `U` or `V`.
  3. You mis-read the row-major layout (all three orders for $x_0$ come first).

---

## Practical tips

| want to see                   | what to do                                                                                                       |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| plain gradient                | set `w = {1,0,0}` and do **not** call `Forward(1, U)`/`Forward(2, V)`                                            |
| Hessian-vector product $H\,U$ | set `w = {0,1,0}` and choose `U`; read the **order 1** block (`dw[3j+1]`)                                        |
| third-tensor contractions     | set `w = {0,0,1}`; `dw[3j+0]` (scaled by 2) gives $f'''(\cdot)$ contracted with `U` and `V` as per formula above |

Always zero-out directions you do **not** want:

```cpp
f.Forward(1 , CppAD::vector<double>(n,0.0));
f.Forward(2 , CppAD::vector<double>(n,0.0));
```

before planting the one non-zero `U` or `V`.

---

### Bottom line

* `U` (first-order direction) and `V` (second-order direction) **become part of the point at which the third-order reverse sweep is evaluated**.
* Even with a diagonal analytic tensor, non-zero `U` can introduce off-diagonal looking numbers in `dw[ * 0 ]` and `dw[ * 1 ]` because those formulas involve products like $f_{jk}U_k$.
* Set `U` and `V` carefully—and interpret the row-major result blocks—to match the derivative you actually need.
