# Final Project 810

## Description

### Human-capital accumulation 
- At any period $t$ agents without a college degree $(\psi=0)$ decide to attend college or not.
- If they attend college, spend $t+1, \ldots, t+5$ periods earning $0$ and after that get a college degree $(\psi = 1)$.

### Firm Side
- Two types of firms $H$ and $L$.
- If employed by $I$ draw wage from $F_I$ with $ \mathbb{E}[H] > \mathbb{E}[L, \psi = 1] > \mathbb{E}[L, \psi = 0]$.
- $H$ firm want to hire above only college graduates.
- $L$ firm want to hire everyone.
    - If $L$ hires you and you are a college graduate then you get a higher wage.

#### Production Function
- Firms have a production function $F_I(h,\psi)$:
$$ f^I(h) = f^I(h, \psi) = A_I( h + \psi\varepsilon), \quad I \in\{L,H\} $$

#### Firm Value Function
$$
\begin{aligned}
J^I_{t}(\omega, \vec{h}) &=(1-\omega) f^I(\vec{h})+\beta_{f} \mathbb{E}\left[(1-\delta) J_{t+1}\left(\omega, \vec{h}^{\prime}\right)\right] \quad \forall t \leq T \\
J^I_{T+1}(\omega, \vec{h}) &=0
\end{aligned}
$$

### Labor market
- Vacancies indexed by $\omega$, $h$ and $I$ -- $v_t(\omega,h,I)$.
- Unemployed workers characterized by $\omega$, $h$, and $\psi$ -- $u_t(\omega,h,\psi)$.
- Labor market tightness $\theta_t(\omega,h,I,\psi)$.
- Job finding rate $p(\cdot) = \frac{M(u,v)}{u}$ and hiring rate $p_f(\cdot) = \frac{M(u,v)}{v}$.
- 
### Agents Value Function

#### Unemployed workers
- Start each period with some level of asset, human capital, and an indicator for having a college degree.
- If they invest in a college degree, their continuation value is the value function 4 years down the road.
- Note that $\psi$ can only move in one direction -- from 0 to 1.
$$
\begin{gathered}
U_{t}(b, \vec{h}, \psi)=\max _{b^{\prime} \geq \underline{b}} u(c)+  \mathbb{E}\left[\max _{\tilde{\omega}, I, \psi'} \beta_{i}(\psi') \left[ p\left(\theta_{t+1}\left(\tilde{\omega}, \vec{h}^{\prime}, I, \psi'\right)\right) W_{t+1}\left(\tilde{\omega}, b^{\prime}, \vec{h}^{\prime}, \psi'\right)+\left(1-p\left(\theta_{t+1}\left(\tilde{\omega}, \vec{h}^{\prime}, I , \psi'\right)\right)\right) U_{t+1}\left(b^{\prime}, \vec{h}^{\prime}, \psi'\right) \right] \right] \forall t \leq T \\
U_{T+1}(b, \vec{h}, \psi)=0
\end{gathered}
$$
where $\beta_{i}(\psi') = \beta^5$ as the continuation value is delayed 5 years into the future.

I am struggling to get a unified notation above, so here is an example.  At each time $t$, if agents do not invest in college then
$$
\begin{gathered}
U_{t}(b, \vec{h}, \psi)=\max _{b^{\prime} \geq \underline{b}} u(c)+ \beta_{i} \mathbb{E}\left[\max _{\tilde{\omega}, I} \left[ p\left(\theta_{t+1}\left(\tilde{\omega}, \vec{h}^{\prime}, I, \psi \right)\right) W_{t+1}\left(\tilde{\omega}, b^{\prime}, \vec{h}^{\prime}, \psi\right)+\left(1-p\left(\theta_{t+1}\left(\tilde{\omega}, \vec{h}^{\prime}, I , \psi\right)\right)\right) U_{t+1}\left(b^{\prime}, \vec{h}^{\prime}, \psi\right) \right] \right] \forall t \leq T \\
U_{T+1}(b, \vec{h}, \psi)=0
\end{gathered}
$$
And, if agents invest in college -- then they have unemployed value function for the next 4 years and then do the labor search in their 5th year.
$$
\begin{gathered}
U_{t}(b, \vec{h}, \psi)=\max _{b^{\prime} \geq \underline{b}} u(c)+ \beta_{i} \mathbb{E}\left[\max _{\tilde{\omega}, I, \psi} \sum_{i=1}^4 U_{t+i} \left(b^{\prime}, \vec{h}^{\prime}, \psi\right) +  \left[ p\left(\theta_{t+5}\left(\tilde{\omega}, \vec{h}^{\prime}, I, \psi \right)\right) W_{t+5}\left(\tilde{\omega}, b^{\prime}, \vec{h}^{\prime}, \psi\right)+\left(1-p\left(\theta_{t+5}\left(\tilde{\omega}, \vec{h}^{\prime}, I , \psi\right)\right)\right) U_{t+5}\left(b^{\prime}, \vec{h}^{\prime}, \psi\right) \right] \right] \forall t \leq T \\
U_{T+1}(b, \vec{h}, \psi)=0
\end{gathered}
$$

We want to think about the tradeoff from having to invest in college  i.e. college degree costs a certain amount for the 4 years -- or just a lump sum when you decide to attend college.

<!-- $$
p(\psi) = p_H(\psi)\mu  + p_L(\psi)(1-\mu)
$$ -->


#### Employed workers
- Might want to quit in search of a college degree. This can make it too complicated!
- Simplest case would be to not allow employed agents to make college decision while employed.

## Simulation
- Agents start $\psi=0$.
- Agents start with different initial wealth.
