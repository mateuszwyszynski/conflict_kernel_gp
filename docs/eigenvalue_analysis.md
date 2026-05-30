# The minimum eigenvalue of a Gaussian kernel on the sphere

This note studies how the minimum eigenvalue of the geodesic Gaussian kernel on $S^2$ behaves,
as a function of the lengthscale $\kappa$ and of the number of sampled points $n$.
Two questions organize the discussion:
how the minimum eigenvalue moves as $\kappa$ ranges from $0$ to $\infty$ (§2),
and how that picture changes as we add more points (§3).
Along the way it becomes negative for a range of $\kappa$,
so the kernel is not positive semidefinite there.

Throughout, $n$ is the number of points,
$\mathbf{1} = (1,\dots,1)^T \in \mathbb{R}^n$ is the all-ones vector,
and $\kappa > 0$ is the lengthscale.

---

## 1. The squared geodesic-distance matrix $D$ and the kernel $K$

Given points $x_1,\dots,x_n \in S^2$, let $d_g(x_i,x_j)$ be the geodesic (great-circle) distance.
The squared geodesic-distance matrix is

$$
D \in \mathbb{R}^{n\times n}, \qquad D_{ij} = d_g(x_i, x_j)^2 .
$$

$D$ is symmetric with zero diagonal ($D_{ii} = 0$) and non-negative entries.
The Gaussian kernel at lengthscale $\kappa$ is the entrywise exponential

$$
K = K(\kappa), \qquad K_{ij} = \exp\!\left(-\frac{D_{ij}}{2\kappa^2}\right).
$$

This is applied element by element; it is not the matrix exponential $e^{-D/2\kappa^2}$.
$K$ is symmetric, has $1$'s on the diagonal, and entries in $(0,1]$.

---

## 2. The behaviour of the minimum eigenvalue as $\kappa$ varies

The quantity of interest is $\lambda_{\min}(K(\kappa))$,
the smallest eigenvalue of the kernel as a function of the lengthscale.
Its overall shape is:

- $\lambda_{\min} \to 1$ as $\kappa \to 0$;
- $\lambda_{\min} \to 0$ as $\kappa \to \infty$;
- $\lambda_{\min} < 0$ starting at some $\kappa^* \in \mathbb{R}$.

The two limits are elementary (§2.1).
The interesting part — that the curve must dip strictly below zero in between —
needs a short perturbative argument (§2.2).

### 2.1 Endpoint behaviour

As $\kappa \to 0$, every off-diagonal entry has $-D_{ij}/(2\kappa^2) \to -\infty$,
so $K_{ij}\to 0$, while the diagonal stays at $K_{ii}=1$. Hence

$$
K(\kappa) \xrightarrow{\ \kappa\to 0\ } I,
$$

whose eigenvalues are all $1$; in particular $\lambda_{\min}\to 1$.

As $\kappa \to \infty$, every exponent $-D_{ij}/(2\kappa^2)\to 0$, so $K_{ij}\to 1$
and

$$
K(\kappa) \xrightarrow{\ \kappa\to\infty\ } J := \mathbf{1}\mathbf{1}^T.
$$

$J = \mathbf{1}\mathbf{1}^T$ has rank one, so $n-1$ of its eigenvalues are $0$;
the single nonzero eigenvalue equals $\operatorname{tr} J = n$, with eigenvector
$\mathbf{1}$. Its zero-eigenspace is $\mathbf{1}^\perp$, the vectors $v$ with
$\mathbf{1}^T v = 0$. So at $\kappa\to\infty$ the minimum eigenvalue of $K$ is
$0$.

### 2.2 The minimum eigenvalue dips below zero in between

We look for some $\kappa$ with $\lambda_{\min}(K) < 0$.

#### 2.2.1 First-order expansion

Write $\epsilon = \dfrac{1}{2\kappa^2}$, so $\epsilon \to 0$ as $\kappa\to\infty$.
Expanding each entry with $e^z = 1 + z + O(z^2)$,

$$
K_{ij} = \exp(-\epsilon D_{ij}) = 1 - \epsilon D_{ij} + O(\epsilon^2),
$$

so $K = J - \epsilon D + O(\epsilon^2)$.

#### 2.2.2 An upper bound on the smallest eigenvalue

We cannot read off $\lambda_{\min}(K)$ by adding the eigenvalues of $J$ and
$-\epsilon D$: in general $\operatorname{eig}(M+N) \neq \operatorname{eig}(M) +
\operatorname{eig}(N)$, because $M$ and $N$ need not share eigenvectors. Instead
we use the Rayleigh-quotient characterization of the smallest eigenvalue,

$$
\lambda_{\min}(K) = \min_{\|v\|=1} v^T K v,
$$

which gives $\lambda_{\min}(K) \le v^T K v$ for every unit vector $v$. It
therefore suffices to exhibit one direction $v^*$ in which the quadratic form is
negative.

#### 2.2.3 Choosing the right vector

We will choose $v^* \in \mathbf{1}^\perp$ with $v^{*T} K v^* < 0$. The orthogonal
projector onto a single vector $w$ is $\dfrac{w w^T}{\|w\|^2}$, so the projector
onto $\operatorname{span}(\mathbf 1)$ is $\frac{1}{n}\mathbf{1}\mathbf{1}^T$ and
the projector onto $\mathbf{1}^\perp$ is

$$
P = I - \frac{1}{n}\mathbf{1}\mathbf{1}^T.
$$

$P$ is symmetric, $P^2 = P$, and $Pv$ is the zero-sum part of $v$ (its entries
sum to $0$); if $v$ is already zero-sum then $Pv = v$.

Consider the symmetric matrix $PDP$, and let $v^*$ be a unit eigenvector for its
largest eigenvalue $\mu := \lambda_{\max}(PDP)$. Since $PDP$ maps into
$\mathbf{1}^\perp$, this $v^*$ is zero-sum, so $P v^* = v^*$ and

$$
v^{*T} D v^* = v^{*T} P D P v^* = \mu .
$$

Evaluating the Rayleigh quotient of $K$ on $v^*$, using $\mathbf 1^T v^* = 0$,

$$
v^{*T} K v^* = \underbrace{(\mathbf 1^T v^*)^2}_{=\,0} - \epsilon\, v^{*T} D v^* + O(\epsilon^2)
= -\,\epsilon\,\mu + O(\epsilon^2),
$$

and combining with the bound above,

$$
\lambda_{\min}(K) \le v^{*T} K v^* = - \epsilon \lambda_{\max}(PDP) + O(\epsilon^2).
$$

Here $D$, and hence $v^*$ and $\mu = \lambda_{\max}(PDP)$, are fixed once the
points are given; they do not depend on $\kappa$. So $v^{*T} K v^*$ is a scalar
function of $\epsilon$ alone with constant Taylor coefficients, $f(\epsilon) =
-\mu\epsilon + c\,\epsilon^2 + \dots$, and the linear term dominates as $\epsilon
\to 0$. Concretely, since $0 \le D_{ij} \le \pi^2$ and $0 \le e^{-x} - (1-x) \le
\tfrac12 x^2$ for $x \ge 0$, the remainder is bounded by $\tfrac12 n \pi^4
\epsilon^2$, so $f(\epsilon) \le -\mu\epsilon + \tfrac12 n \pi^4 \epsilon^2$,
which is negative whenever $\epsilon < 2\mu / (n\pi^4)$.

So if $\mu > 0$, then for sufficiently small $\epsilon$ (large $\kappa$) we have
$\lambda_{\min}(K) < 0$. Together with $\lambda_{\min} = 1$ near $\kappa = 0$ and
$\lambda_{\min}\to 0$ as $\kappa\to\infty$, the minimum eigenvalue must cross
zero and spend an interval of $\kappa$ negative — the dip seen in the plot — then
recover toward $0$, since the bound $-\epsilon\mu \to 0^-$.

---

## 3. How the crossing depends on the number of points $n$

Section 2 showed that the dip below zero is driven by the existence of a
zero-sum direction $v$ with $v^T D v > 0$. Whether such a direction exists, and
where the crossing lands, depends on how many points we sample. Write
$K^{(n)}(\kappa)$ for the kernel of an $n$-point configuration and define the
crossing point as the onset of negativity,

$$
\kappa^*_n = \inf\{\,\kappa > 0 : \lambda_{\min}(K^{(n)}(\kappa)) < 0\,\},
$$

with $\kappa^*_n = +\infty$ when the kernel stays positive for every $\kappa$
(no crossing). Two facts shape the picture: a crossing is impossible until
$n = 4$ (§3.1), and once it exists it only moves to smaller $\kappa$ as points
are added (§3.2).

### 3.1 A crossing requires at least $n = 4$ points

For $n \le 3$ there is no crossing: $\lambda_{\min}(K^{(n)}(\kappa)) \ge 0$ for
every $\kappa$.

The reason is that three or fewer points can always be laid out flat without
distorting their distances. Two points lie on a line; three points with pairwise
distances $d_{12}, d_{13}, d_{23}$ obey the triangle inequality (a geodesic
distance is a genuine metric), so they form an ordinary triangle in the plane
with those side lengths. We can therefore place points $p_1,\dots,p_n \in
\mathbb{R}^2$ with $D_{ij} = \|p_i - p_j\|^2$. But then $K$ is just the familiar
Gaussian kernel on flat points,

$$
K_{ij} = \exp\!\left(-\frac{\|p_i - p_j\|^2}{2\kappa^2}\right),
$$

which is positive definite for every bandwidth $\kappa$. So the kernel never goes
negative, and $\kappa^*_n = +\infty$. Consistently, for any zero-sum $v$,

$$
v^T D v = \sum_{i,j} v_i v_j\big(\|p_i\|^2 + \|p_j\|^2 - 2\,p_i^T p_j\big)
        = -2\,\Big\|\sum_i v_i p_i\Big\|^2 \le 0,
$$

so the direction needed in §2 simply does not exist.

This breaks at $n = 4$. Four points carry six pairwise distances, and a curved
sphere generally cannot reproduce all six with any flat configuration — the
distances no longer come from Euclidean points, and the kernel can lose
positivity. It does not happen for every configuration (empirically only about
half of random $4$-point samples on $S^2$ cross), but $n = 4$ is the smallest
size at which it can.

### 3.2 The crossing point decreases as $n$ grows

Adding a point appends one row and column to the kernel, so $K^{(n)}(\kappa)$ is
a principal submatrix of $K^{(n+1)}(\kappa)$. The smallest eigenvalue of a
symmetric matrix is never larger than that of any principal submatrix (Cauchy
interlacing), so

$$
\lambda_{\min}(K^{(n+1)}(\kappa)) \le \lambda_{\min}(K^{(n)}(\kappa))
\qquad \text{for every } \kappa .
$$

At fixed $\kappa$, adding points can only push $\lambda_{\min}$ down. In
particular, if the $n$-point kernel is already negative at some $\kappa$, so is
the $(n+1)$-point kernel — the set of $\kappa$ where the kernel is negative only
grows. Taking the infimum of that set,

$$
\kappa^*_{n+1} \le \kappa^*_n .
$$

For nested configurations the crossing point is monotonically non-increasing in
$n$. The same holds on average: an $n$-point sample has the same distribution as
an $(n+1)$-point sample with one point removed, and removing a point only raises
the threshold, so $\mathbb{E}[\kappa^*_{n+1}] \le \mathbb{E}[\kappa^*_n]$.

### 3.3 The limit of many points

As $n \to \infty$ the points fill the sphere and the kernel matrix approaches a
smoothing operator that averages a function against $\exp(-d_g(x,y)^2/2\kappa^2)$.
Because the kernel depends only on the angle between $x$ and $y$, this operator
is diagonalized by spherical harmonics, and its eigenvalues can be computed one
frequency at a time. One finds a negative eigenvalue for every $\kappa > 0$ (at
small $\kappa$ it lives at high frequency and is extremely small). So in the
limit the geodesic Gaussian is never positive for any bandwidth, and there is no
positive limiting threshold: $\mathbb{E}[\kappa^*_n] \to 0$.
