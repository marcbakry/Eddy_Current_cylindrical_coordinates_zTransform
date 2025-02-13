## Courants de Foucault transient en coordonnées cylindriques avec symétrie de rotation autour de l'axe $\vec e_z$

On a l'équation en potentiel $\mathbf A$

$$
\nabla\times\left(\nu\nabla\times\mathbf A\right) + \sigma\dfrac{\partial\mathbf A}{\partial t} = \mathbf J
$$

puis

$$\mathbf B = \nabla\times\mathbf A.$$

En coordonnées cylindriques $(\rho,\theta,z)$, l'opérateur rotationnel est tel que 

$$
\nabla\times\mathbf A = \begin{pmatrix}
    \dfrac{1}{\rho}\dfrac{\partial A_z}{\partial\theta} - \dfrac{\partial A_{\theta}}{\partial z} \\
    \dfrac{\partial A_{\rho}}{\partial z} - \dfrac{\partial A_z}{\partial\rho} \\
    \dfrac{1}{\rho}\dfrac{\partial (\rho A_{\theta})}{\partial\rho} - \dfrac{1}{\rho}\dfrac{\partial A_{\rho}}{\partial\theta}
\end{pmatrix}
$$

Par symétrie, $\dfrac{\partial}{\partial\theta}\equiv 0$ donc

$$
\nabla\times\mathbf A = \begin{pmatrix}
-\dfrac{\partial A_{\theta}}{\partial z} \\
\dfrac{\partial A_{\rho}}{\partial z} - \dfrac{\partial A_z}{\partial\rho} \\
\dfrac{1}{\rho}\dfrac{\partial (\rho A_{\theta})}{\partial\rho}
\end{pmatrix}
$$

Par ailleurs, $\mathbf J = J_{\theta}\,\vec e_{\theta}$ donc $\mathbf A = A_{\theta}\,\vec e_{\theta}$ ce qui implique que

$$
\nu\nabla\times\mathbf A = \nu\begin{pmatrix}
-\dfrac{\partial A_{\theta}}{\partial z} \\ 
0 \\
\dfrac{1}{\rho}\dfrac{\partial(\rho A_{\theta})}{\partial\rho}
\end{pmatrix}.
$$
On va supposer que $\nu$ est constant par morceaux. Par conséquent,

$$
\nabla\times\left(\nu\nabla\times\mathbf A\right)=\begin{pmatrix}
0 \\
-\dfrac{\partial^2 A_{\theta}}{\partial z^2} - \dfrac{\partial}{\partial\rho}\left(\dfrac{1}{\rho}\dfrac{\partial (\rho A_{\theta})}{\partial\rho}\right)\\
0
\end{pmatrix}.
$$

<!-- Par ailleurs, $\partial_{\rho}(\rho A_{\theta}) = \rho\partial_{\rho}A_{\theta} + A_{\theta}$ donc

$$\dfrac{1}{\rho}\partial_{\rho}(\rho A_{\theta}) = \partial_{\rho}A_{\theta}+\dfrac{A_{\theta}}{\rho}$$

puis

$$
\partial_{\rho}{\left(\dfrac{1}{\rho}\partial_{\rho}(\rho A_{\theta})\right)} = \partial_{\rho\rho}A_{\theta} + \dfrac{1}{\rho}\partial_{\rho}A_{\theta} - \dfrac{1}{\rho^2}A_{\theta}.
$$

Finalement, on est amené à résoudre

$$
\nu\left(-\partial_{\rho\rho}A_{\theta}-\partial_{zz}A_{\theta}-\dfrac{1}{\rho}\partial_{\rho}A_{\theta}+\dfrac{1}{\rho^2}A_{\theta}\right) + \sigma \partial_t A_{\theta} = J_{\theta}
$$ -->

On doit donc résoudre

$$
\nu\left(-\dfrac{\partial^2 A_{\theta}}{\partial z^2} - \dfrac{\partial}{\partial\rho}\left(\dfrac{1}{\rho}\dfrac{\partial (\rho A_{\theta})}{\partial\rho}\right)\right) + \sigma \partial_t A_{\theta} = J_{\theta}
$$

Testons cette équation par une fonction $A'_{\theta}$ et intégrons sur le domaine $\Omega$ en supposant une condition de Dirichlet homogène sur le bord du domaine de calcul, il vient

$$
\begin{aligned}
\int_{\Omega}{\left(-A'_{\theta}\partial_{zz}A_{\theta} - A'_{\theta}\partial_{\rho}\left(\dfrac{1}{\rho}\partial_{\rho}\left(\rho A_{\theta}\right)\right)\right)\,\rho d\rho d\theta dz} &= \int_{\Omega}{\left(\partial_z A'_{\theta}\partial_{z}A_{\theta} - A'_{\theta}\partial_{\rho}\left(\dfrac{1}{\rho}\partial_{\rho}\left(\rho A_{\theta}\right)\right)\right)\,\rho d\rho d\theta dz}, \\
&= \int_{\Omega}{\left(\partial_z A'_{\theta}\partial_{z}A_{\theta} + \dfrac{1}{\rho^2} \partial_{\rho}\left(\rho A'_{\theta}\right)\partial_{\rho}\left(\rho A_{\theta}\right)\right)\,\rho d\rho d\theta dz}
\end{aligned}
$$

Ensuite,
$$\partial_\rho (\rho A_{\theta}) = A_{\theta} + \rho\partial_{\rho}A_{\theta}$$

Après résolution, le champ $\mathbf B$ sur l'axe $\rho=0$ peut être calculé par extrapolation en calculant le champ en $\rho=\varepsilon\ll 1$.

#### Résolution d'un problème temporel par $z$-transform

L'équation transiente en haut de ce document devient, après application de la $z$-transform, 

$$\nabla\times\left(\nu\nabla\times\mathbf A_z\right) + \dfrac{\sigma(1 - z)}{\Delta t}\mathbf A_z = \mathbf J_z$$

qui est l'analogue discret d'un problème harmonique de fréquence

$$\omega = \dfrac{\sigma (1-z)}{i\Delta t}$$