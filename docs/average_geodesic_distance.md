# Average geodesic distance between two random points on the sphere

Pick two points independently and uniformly on the unit sphere $S^2$. By
rotational symmetry, fix the first point at the north pole; the geodesic
(great-circle) distance to the second is then its polar angle $\theta \in [0,
\pi]$. The average distance is therefore $\mathbb{E}[\theta]$, and we need the
distribution of $\theta$.

"Uniform on the sphere" means uniform with respect to surface area, not uniform
in the angle $\theta$. With polar angle $\theta \in [0,\pi]$ and azimuth $\varphi
\in [0, 2\pi)$, the area element is

$$
dA = \sin\theta \, d\theta \, d\varphi .
$$

The $\sin\theta$ factor reflects that a band at polar angle $\theta$ is a circle
of radius $\sin\theta$, so there is more area near the equator than near the
poles. The total area is $4\pi$, and the uniform-on-area joint density is

$$
p(\theta, \varphi) = \frac{\sin\theta}{4\pi}, \qquad \theta\in[0,\pi],\ \varphi\in[0,2\pi).
$$

Integrating out $\varphi$ gives the marginal density of $\theta$,

$$
p(\theta) = \int_0^{2\pi} \frac{\sin\theta}{4\pi}\, d\varphi = \frac{\sin\theta}{2}, \qquad \theta \in [0,\pi],
$$

which integrates to $1$ and is symmetric about $\theta = \pi/2$. The mean is

$$
\mathbb{E}[\theta] = \int_0^\pi \theta \cdot \tfrac{1}{2}\sin\theta \, d\theta = \frac{\pi}{2},
$$

using $\int_0^\pi \theta\sin\theta\,d\theta = \big[-\theta\cos\theta\big]_0^\pi +
\int_0^\pi \cos\theta\,d\theta = \pi$. The symmetry of $p(\theta)$ about $\pi/2$
gives the same answer directly.

## A note on sampling uniformly

Area-uniform sampling requires drawing $\theta$ with density $\tfrac12\sin\theta$,
e.g. $\theta = \arccos(1 - 2u)$ with $u \sim \mathrm{Unif}[0,1]$. Drawing
$\theta \sim \mathrm{Unif}[0,\pi]$ instead is uniform in the angle, not on the
sphere, and over-samples the poles.
