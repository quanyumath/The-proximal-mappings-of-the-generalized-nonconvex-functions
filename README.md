# The proximal mappings of the generalized nonconvex functions

## Generalized nonconvex functions
We list five of them as follows. Here we only consider the case $x > 0$.
- L1: $\psi^{\text{L1}}\left(x\right)=x$;
- Lp: $\psi^{\text{Lp}}\left(x\right)=x^p,~p\in \left(0,1\right)$;
- MCP: 
```math
\psi^{\text{MCP}}(x)=\begin{cases}x-\frac{x^2}{2\alpha}, & 0 \leq x \leq \alpha, \\ \frac{\alpha}{2}, & x>\alpha \end{cases}
```
- Logarithm: $\psi^{\text{Log}}(x)=\log (\frac{x}{\theta}+1)$ with $\theta>0$;
- Capped folded functions:
  - Capped $L1$: $\psi^{\text{CapL1}}\left(x\right)=\min\left\lbrace 1,\frac{x}{v}\right\rbrace$; 
  - Capped $Lp$: $\psi^{\text{CapLp}}\left(x\right)=\min\left\lbrace 1,\frac{x^p}{v^p}\right\rbrace ,~p\in \left(0,1\right)$; 
  - Capped MCP: $\psi^{\text{CapMCP}}\left(x\right)=\min\left\lbrace1, \frac{2\alpha}{\nu(2\alpha-\nu)} \psi^{\text{MCP}}(x)\right\rbrace, ~0<\nu<\alpha$;
  - Capped Logarithm: $\psi^{\text{CapLog}}\left(x\right)=\min\left\lbrace1,\frac{1}{\psi^{\text{Log}}(v)}\psi^{\text{Log}}(x)\right\rbrace.$


注意：prox_tnn_psi.m 和 Psi.m 中 \psi 的参数 nv, p, alpha, theta 需要保持一致

