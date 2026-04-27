To compute the Hessian applied to a perturbation $\Delta p$ above, we need to compute the FP, AP, second-order FP and second-order AP. 

Isolating the $p$-variation, we get:

$$
\begin{aligned}
    \partial_p \mathcal{L}_2 &= 
    - \left( V_j^\dagger \frac{\partial^2 F_j}{\partial p^2} \varepsilon_p + \Psi^\dagger \frac{\partial^2 G}{\partial p^2} \varepsilon_p, \delta p \right)_{\Gamma \times \mathcal{T}} \\
    &\quad - \left( \Delta V_j^\dagger \frac{\partial F_j}{\partial p}, \delta p \right)_{\Gamma \times \mathcal{T}} - \left( \Delta \Psi^\dagger \frac{\partial G}{\partial p}, \delta p \right)_{\Gamma \times \mathcal{T}} \\
    &\quad - \left( V_j^\dagger \left[ \frac{\partial^2 F_j}{\partial V_i \partial p} \Delta V_i + \frac{\partial^2 F_j}{\partial \Psi \partial p} \Delta \Psi \right], \delta p \right)_{\Gamma \times \mathcal{T}} \\
    &\quad - \left( \Psi^\dagger \left[ \frac{\partial^2 G}{\partial V_i \partial p} \Delta V_i + \frac{\partial^2 G}{\partial \Psi \partial p} \Delta \Psi \right], \delta p \right)_{\Gamma \times \mathcal{T}}
\end{aligned}
$$

---

### Forward Problem

$$
\begin{aligned}
    \rho\ddot{u}_{J} &= \partial_{I}C_{IJKL}\partial_{K}u_{L} + Q_{J}, && \overline{x}\in\Omega, && t\in\mathcal{T}, \\
    \overline{u} &= \overline{u}_{0}, \quad \dot{\overline{u}} = \overline{v}_{0}, && \overline{x}\in\Omega, && t=0, \\
    L\overline{u} &= \overline{g}, && \overline{x}\in\partial\Omega, && t\in\mathcal{T}, \\
    V_{n} &= 0, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \overline{T}^{+} &= -\overline{T}^{-}, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \tau_{J} &= F_{J}(\overline{V}, \Psi, p), && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \dot{\Psi} &= G(\overline{V}, \Psi, p), && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \Psi &= \Psi_{0}, && \overline{x}\in\Gamma, && t=0.
\end{aligned}
$$

---

### Adjoint Problem

$$
\begin{aligned}
    \rho\ddot{u}_{J}^{\dagger} &= \partial_{I}C_{IJKL}\partial_{K}u_{L}^{\dagger} + Q_{J}^{\dagger}, && \overline{x}\in\Omega, && t\in\mathcal{T}, \\
    \overline{u}^{\dagger} &= \overline{0}, \quad \dot{\overline{u}}^{\dagger} = \overline{0}, && \overline{x}\in\Omega, && t=T, \\
    L^{\dagger}\overline{u}^{\dagger} &= 0, && \overline{x}\in\partial\Omega, && t\in\mathcal{T}, \\
    V_{n}^{\dagger} &= 0, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \overline{T}^{\dagger+} &= -\overline{T}^{\dagger-}, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \tau_{J}^{\dagger} &= F_{J}^{\dagger}(\overline{V}^{\dagger}, \Psi^{\dagger}), && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    -\dot{\Psi}^{\dagger} &= G^{\dagger}(\overline{V}^{\dagger}, \Psi^{\dagger}), && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \Psi^{\dagger} &= 0, && \overline{x}\in\Gamma, && t=T.
\end{aligned}
$$

---

### Second-Order Forward Problem

$$
\begin{aligned}
    \rho\Delta\ddot{u}_{J} &= \partial_{I}C_{IJKL}\partial_{K}\Delta u_{L}, && \overline{x}\in\Omega, && t\in\mathcal{T}, \\
    \Delta\overline{u} &= \overline{0}, \quad \Delta\dot{\overline{u}} = \overline{0}, && \overline{x}\in\Omega, && t=0, \\
    L\Delta\overline{u} &= \overline{0}, && \overline{x}\in\partial\Omega, && t\in\mathcal{T}, \\
    \Delta V_{n} &= 0, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \Delta\overline{T}^{+} &= -\Delta\overline{T}^{-}, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \Delta\tau_{J} &= \Delta F_{J}, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \Delta\dot{\Psi} &= \Delta G, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \Delta\Psi &= 0, && \overline{x}\in\Gamma, && t=0.
\end{aligned}
$$

$$
\begin{aligned}
    \Delta F_{j} &= \frac{\partial G^\dagger}{\partial V_i^\dagger} \Delta \Psi + \frac{\partial F_j^\dagger}{\partial V_i^\dagger} \Delta V_i + \frac{\partial F_j}{\partial p}\varepsilon_p \\
    \Delta G &= -\frac{\partial G^\dagger}{\partial \Psi^\dagger}\Delta \Psi - \frac{\partial F_j^\dagger}{\partial \Psi^\dagger}\Delta V_j + \frac{\partial G}{\partial p}\varepsilon_p
\end{aligned}
$$

---

### Second-Order Adjoint Problem

$$
\begin{aligned}
    \rho\Delta\ddot{u}_{J}^{\dagger} &= \partial_{I}C_{IJKL}\partial_{K}\Delta u_{L}^{\dagger} + \Delta Q_{J}^{\dagger}, && \overline{x}\in\Omega, && t\in\mathcal{T}, \\
    \Delta\overline{u}^{\dagger} &= \overline{0}, \quad \Delta\dot{\overline{u}}^{\dagger} = \overline{0}, && \overline{x}\in\Omega, && t=T, \\
    L^{\dagger}\Delta\overline{u}^{\dagger} &= 0, && \overline{x}\in\partial\Omega, && t\in\mathcal{T}, \\
    \Delta V_{n}^{\dagger} &= 0, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \Delta\overline{T}^{\dagger+} &= -\Delta\overline{T}^{\dagger-}, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \Delta\tau_{J}^{\dagger} &= \Delta F_{J}^{\dagger}, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    -\Delta\dot{\Psi}^{\dagger} &= \Delta G^{\dagger}, && \overline{x}\in\Gamma, && t\in\mathcal{T}, \\
    \Delta\Psi^{\dagger} &= 0, && \overline{x}\in\Gamma, && t=T.
\end{aligned}
$$

$$
\begin{aligned}
    \Delta F_{J}^{\dagger} &= \frac{\partial G^\dagger}{\partial V_i} \Delta \Psi + \frac{\partial F_j^\dagger}{\partial V_i} \Delta V_i + \frac{\partial G}{\partial V_i} \Delta \Psi^\dagger + \frac{\partial F_j}{\partial V_i} \Delta V_i^\dagger \\
    &\quad + \varepsilon_p \left( \Psi^\dagger \frac{\partial^2 G}{\partial V_i \partial p} + V_j^\dagger\frac{\partial^2 F_j}{\partial V_i \partial p} \right), \\
    \Delta G^{\dagger} &= -\frac{\partial G^\dagger}{\partial \Psi}\Delta \Psi - \frac{\partial F_j^\dagger}{\partial \Psi}\Delta V_j -\frac{\partial G}{\partial \Psi}\Delta \Psi^\dagger - \frac{\partial F_j}{\partial \Psi}\Delta V_j^\dagger \\
    &\quad + \varepsilon_p \left( \Psi^\dagger \frac{\partial^2 G}{\partial \Psi \partial p} + V_j^\dagger\frac{\partial^2 F_j}{\partial \Psi \partial p} \right)
\end{aligned}
$$

$$
    \Delta Q^\dagger_j = \sum_{k=1}^{N_{rec}} \Delta \dot{u}_j \hat{\delta}(\overline{x} - \overline{x}_r^{(k)})
$$

So the second-order forward problem drives the second-order adjoint problem with $\Delta Q_j^\dagger$.