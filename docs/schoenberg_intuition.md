# Intuition for Schoenberg's theorem

This companion note builds intuition for *why* conditional negative definiteness
(CND) of the squared distance is exactly the condition that makes a Gaussian
kernel positive semidefinite at **all** lengthscales, and why the sphere's
geodesic metric fails it.

---

## 0. Statement

> **Schoenberg's theorem (1938).** For a symmetric matrix $D$ with zero diagonal
> (a "squared-distance matrix"), the Gaussian kernel
> $K_{ij} = \exp(-t\,D_{ij})$ is **positive semidefinite for every $t > 0$** if
> and only if $D$ is **conditionally negative definite** — equivalently, if and
> only if the distances embed isometrically into a Hilbert space, i.e. there are
> vectors $\varphi(x_i)$ in some Hilbert space with
> $D_{ij} = \lVert \varphi(x_i) - \varphi(x_j)\rVert^2$.

Here $t = 1/(2\kappa^2)$ in the lengthscale language of the other note.

---

## 1. "PSD kernel" — yes, just positive semidefinite

Yes. A kernel matrix $K$ is **positive semidefinite (PSD)** if it is symmetric
and

$$
v^\top K v \ge 0 \qquad \text{for all } v \in \mathbb{R}^n,
$$

equivalently all its eigenvalues are $\ge 0$. "Positive *definite*" would be the
strict version ($> 0$ for $v\neq 0$). In machine learning, PSD is the property
that makes $K$ a valid Gram matrix of inner products — i.e. a legitimate kernel.
A **negative** eigenvalue (what the plot shows) means $K$ is *not* PSD and hence
*not* a valid kernel for that lengthscale.

---

## 2. What "a Hilbert space" means here

It is a **general** (real) Hilbert space — not necessarily $\mathbb{R}^d$.

> **Reminder — Hilbert space.** A Hilbert space $H$ is a vector space equipped
> with an inner product $\langle\cdot,\cdot\rangle$ (giving a notion of angle and
> length, $\lVert u\rVert = \sqrt{\langle u,u\rangle}$) that is **complete**:
> every Cauchy sequence converges in $H$. Finite-dimensional examples are just
> $\mathbb{R}^d$ with the dot product. Infinite-dimensional examples include
> $\ell^2$ (square-summable sequences) and $L^2$ (square-integrable functions).
> Completeness is the technical condition that lets you do analysis (limits,
> projections, orthonormal bases) exactly as in $\mathbb{R}^d$.

For a *finite* set of $n$ points you never need more than $n$ dimensions, so
"embeds in a Hilbert space" effectively means "embeds in some $\mathbb{R}^m$ with
the ordinary Euclidean inner product." The Hilbert-space language is just the
clean, dimension-free way to state it.

---

## 3. Why CND is exactly the right condition — the intuition

Two complementary ways to see it.

### 3a. PSD kernel ⟺ feature map ⟺ Euclidean-like squared distance

$K$ is PSD if and only if there is a feature map $\varphi$ with
$K_{ij} = \langle \varphi(x_i), \varphi(x_j)\rangle$. For a Gaussian kernel,
$K_{ij} = \exp(-t\,D_{ij})$, so asking "$K$ PSD for all $t$" is asking that
$D_{ij}$ behave like a **squared Euclidean distance** in some feature space:
$D_{ij} = \lVert\varphi(x_i) - \varphi(x_j)\rVert^2$. The set of squared
distances realizable by *some* point configuration in a Hilbert space is
characterized **exactly** by conditional negative definiteness. So:

$$
\text{geometry embeds in a Hilbert space} \iff D \text{ is CND} \iff \text{Gaussian kernel PSD for all } t .
$$

CND is "flat enough to be laid down in a flat (Euclidean) space."

### 3b. The small-$t$ limit is literally the argument from the other note

This is the most concrete intuition, and it ties the two notes together. Expand
for small $t$ (= large lengthscale):

$$
K = \mathbf 1\mathbf 1^\top - t\,D + O(t^2).
$$

On the zero-sum subspace $\mathbf 1^\perp$, the dominant term $\mathbf 1\mathbf
1^\top$ contributes nothing, so PSD-ness there reduces, as $t\to 0$, to

$$
v^\top K v \approx -t\, v^\top D v \ge 0 \quad\text{for all zero-sum } v
\iff v^\top D v \le 0 \quad\text{for all zero-sum } v,
$$

which is **precisely the definition of $D$ being CND**. So the small-$t$ edge of
Schoenberg's theorem is exactly the perturbation argument in
`eigenvalue_analysis.md` §2.2. Schoenberg's real content is the remarkable fact
that this single condition, read off at the small-$t$ edge, actually guarantees
PSD-ness for **every** $t>0$ at once — not just in the limit. If $A$ is *not*
CND, the same edge already exhibits a negative direction, and PSD-ness fails for
a whole range of $t$.

---

## 4. Why the sphere's geodesic metric fails — and the cartography connection

The sphere with **geodesic** distance does not embed isometrically into any
Hilbert space: there is no assignment of vectors $\varphi(x)$ in a flat space
whose Euclidean distances reproduce great-circle distances. Intuitively, the
sphere is **intrinsically curved** (positive Gaussian curvature), and a flat
space has zero curvature; you cannot match the two without distortion.

> **Yes — this is the same circle of ideas as "no perfect flat map of the
> Earth."** The cartographic fact you are remembering is a consequence of
> **Gauss's *Theorema Egregium***: Gaussian curvature is *intrinsic* — preserved
> by any distance-preserving (isometric) map. The sphere has curvature
> $1/R^2 > 0$; a flat plane (or a flat Hilbert space) has curvature $0$. Since an
> isometry must preserve curvature, **no isometric map sphere → plane exists**.
> That is exactly why every world map distorts *something*: Mercator preserves
> angles (conformal) but blows up areas near the poles; equal-area projections
> preserve areas but distort angles/shapes; and **no** projection preserves the
> actual geodesic *distances*. "Geodesic distances can't be laid flat without
> distortion" is the same statement as "the geodesic squared-distance matrix is
> not CND / doesn't embed in a Hilbert space."

A small concrete witness of non-embeddability: take three points and form a thin
spherical triangle near a pole, or use near-antipodal points (geodesic distance
up to $\pi$). The curvature makes some triple of geodesic distances violate the
relations any flat (Euclidean) triangle must satisfy. Equivalently, there is a
zero-sum vector $v$ with $v^\top D v > 0$ — a direction certifying $D$ is not
CND, which is the algebraic shadow of the geometric obstruction.

By contrast, the **chordal** (straight-line, through-the-ball) distance on the
sphere is, by construction, the Euclidean distance of points already sitting in
$\mathbb{R}^3$. Its squared-distance matrix *is* CND, so the chordal Gaussian
kernel stays PSD for all lengthscales. The negativity in the plot is therefore a
specific consequence of insisting on the **intrinsic geodesic** metric.

---

## 5. One-line summary

A Gaussian kernel is a valid (PSD) kernel at all scales exactly when its squared
distance is "flat" (CND / Hilbert-embeddable). The sphere's geodesic distance is
intrinsically curved — the same reason no flat map of the Earth can preserve
distances — so by Schoenberg's theorem its Gaussian kernel must develop negative
eigenvalues for some lengthscales.
