# Final Project 810

## Description

### Human-capital accumulation 
- At any period $t$ agents without a college degree $(\psi=0)$ decide to attend college or not.
- If they attend college, spend $t+1, \ldots, t+5$ periods earning $0$ and after that get a college degree $(\psi = 1)$.


#### Agents Value Function

$$
\begin{gathered}
U_{t}(b, \vec{h}, \psi)=\max _{b^{\prime} \geq \underline{b}} u(c)+\beta_{i} \mathbb{E}\left[\max _{\tilde{\omega}} p\left(\theta_{t+1}\left(\tilde{\omega}, \vec{h}^{\prime}\right)\right) W_{t+1}\left(\tilde{\omega}, b^{\prime}, \vec{h}^{\prime}\right)+\left(1-p\left(\theta_{t+1}\left(\tilde{\omega}, \vec{h}^{\prime}\right)\right)\right) U_{t+1}\left(b^{\prime}, \vec{h}^{\prime}\right)\right] \forall t \leq T \\
U_{T+1}(b, \vec{h})=0
\end{gathered}
$$

$$
p(\psi) = p_H(\psi)\mu  + p_L(\psi)(1-\mu)
$$
### Firm Side
- Two types of firms $H$ and $L$.
- If employed by $I$ draw wage from $F_I$ with $\mathbb{E}[H]>\mathbb{E}[L, \psi = 1] >\mathbb{E}[L, \psi = 0] $.
- $H$ firm want to hire above only college graduates.
- $L$ firm want to hire everyone.
    - If $L$ hires you and you are a college graduate then you get a higher wage.

#### Production Function
- Firms have a production function $F_I(h,\psi)$:
$$ F_I(h, \psi) = A_I( h + \psi\varepsilon), \quad I \in\{L,H\} $$

#### Firm Value Function
$$
\begin{aligned}
J^H_{t}(\omega, \vec{h}) &=(1-\omega) f(\vec{h})+\beta_{f} \mathbb{E}\left[(1-\delta) J_{t+1}\left(\omega, \vec{h}^{\prime}\right)\right] \quad \forall t \leq T \\
J^H_{T+1}(\omega, \vec{h}) &=0
\end{aligned}
$$


## Simulation
- Agents start $z=0$.
- Agents start with different initial wealth.
