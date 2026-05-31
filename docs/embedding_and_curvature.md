# Embedding point sets in Euclidean space, and why curvature obstructs it

In §3 of `eigenvalue_analysis.md` the crossing was tied to whether the squared
geodesic-distance matrix $D$ comes from a Euclidean configuration: if we can
place points $p_1,\dots,p_n$ in some Euclidean space with $D_{ij} = \|p_i -
p_j\|^2$, the kernel stays positive for every $\kappa$; otherwise it eventually
goes negative. This note looks at that embedding question geometrically — what
space the points would live in, and why a curved sphere keeps the geodesic
distances from fitting into a flat one.

## 1. Why $n$ points embed in $\mathbb{R}^{n-1}$

Suppose $D$ does come from Euclidean points $p_1,\dots,p_n$. Translating so that
$p_1$ sits at the origin, the configuration is spanned by the $n-1$ vectors
$p_2 - p_1, \dots, p_n - p_1$, so it lies in a subspace of dimension at most
$n-1$. That is the ceiling, and each new point can raise the dimension by at most
one: three points span a plane, a fourth reaches into a third dimension, and so
on. There is no fixed low-dimensional space the points must fit into — the room
grows with $n$.

Whether such points exist at all is decided by a single matrix. Center $D$ with
the projector $P = I - \frac1n \mathbf 1\mathbf 1^T$ from §2 and form

$$
G = -\tfrac12\, P D P .
$$

If $D_{ij} = \|p_i - p_j\|^2$, this $G$ is exactly the Gram matrix $G_{ij} =
(p_i - \bar p)^T (p_j - \bar p)$ of the centered points. Conversely, the points
exist precisely when $G$ is positive semidefinite: writing $G = X X^T$ recovers
coordinates $X \in \mathbb{R}^{n \times r}$ with $r = \operatorname{rank}(G)$,
and $r \le n-1$ because $G\mathbf 1 = 0$. So "embeddable" means "$G \succeq 0$",
and the embedding dimension is $\operatorname{rank}(G)$.

This is the §3 condition seen through $G$ instead of $PDP$. Since $PDP$ always
has $\mathbf 1$ in its kernel, $G \succeq 0$ is the same as $\lambda_{\max}(PDP)
= 0$, and a positive $\lambda_{\max}(PDP)$ is a negative eigenvalue of $G$ — a
configuration no Euclidean point set can produce.

## 2. Why curvature obstructs embedding

On the sphere the distance between two points is the length of the great-circle
arc joining them — a curved path. Arc length is always at least the straight-line
(chord) distance, and the gap widens as the points separate. To flatten the
configuration we would need a Euclidean space in which these arc lengths are
straight-line distances, but straight lines are exactly what the arcs are not.
This is Gauss's Theorema Egregium: curvature is intrinsic to a surface, so a
distance-preserving map from a curved surface into a flat space cannot exist —
the same reason no flat map of the Earth preserves all distances.

A concrete witness. Take four points on the equator, $90^\circ$ apart. Their arc
distances are $\pi/2$ between neighbours and $\pi$ between opposite points. A flat
realization would need a quadrilateral with all four sides $\pi/2$; placed as a
square, its diagonal is only $\frac{\pi}{2}\sqrt2 \approx 2.22$, whereas the
sphere demands $\pi \approx 3.14$. The arc that wraps around the sphere is longer
than flat geometry allows for those side lengths, so no Euclidean placement
reproduces all six distances and $G$ acquires a negative eigenvalue.

Two features of the picture follow from the same idea.

It is mild for nearby points and only appears once $n \ge 4$. Over a small cap
the sphere is nearly flat — arc and chord almost agree — so a few close points
embed with only a tiny error, and three points always embed exactly (any triangle
is flat). The discrepancy grows with separation and with the number of points:
embeddability requires every four-point subset to flatten at once, and
arrangements strung along great circles, like the equatorial example, are the
ones that fail hardest.

It is not raw spread that breaks embedding. Four points at the vertices of a
regular tetrahedron are maximally spread yet embed perfectly, because any set of
equidistant points realizes a regular simplex. What matters is whether the
arrangement bends around the sphere's curved geodesics, not how far apart the
points are.
