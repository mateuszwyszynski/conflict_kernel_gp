"""
Streamlit app: minimum eigenvalue of the geodesic Gaussian (RBF) kernel on S^2.

Interactively control the number of points per sample, the number of samples,
the lengthscale grid, and the sphere-sampling scheme. The app reproduces the
"minimum eigenvalue vs lengthscale" figure and additionally computes the
**average zero-crossing lengthscale** (where the minimum eigenvalue first turns
negative), per sample and for the mean curve.

Run with:
    uv run streamlit run streamlit_app.py
"""

import math

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st


# --------------------------------------------------------------------------- #
# Geometry / kernel helpers
# --------------------------------------------------------------------------- #
def geodesic_arc(u, v):
    """Great-circle distance (arc-length) on S^2 between unit vectors u, v."""
    dot = float(np.dot(u, v))
    dot = max(-1.0, min(1.0, dot))  # numerical safety
    return math.acos(dot)


def build_squared_distance_matrix(points_xyz):
    """Pairwise squared geodesic-distance matrix D_ij = d_g(x_i, x_j)^2."""
    n = len(points_xyz)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d2 = geodesic_arc(points_xyz[i], points_xyz[j]) ** 2
            D[i, j] = d2
            D[j, i] = d2
    return D


def gaussian_gram(D, kappa):
    """Gaussian Gram from squared geodesic-distance matrix D (entrywise exp)."""
    return np.exp(-D / (2 * kappa**2))


def generate_sphere_points(n_points, n_samples, rng, uniform_on_sphere=True):
    """Generate `n_samples` collections of `n_points` points on S^2.

    uniform_on_sphere=True  -> correct area-uniform sampling: theta=arccos(1-2u)
                               so that cos(theta) ~ Unif[-1, 1].
    uniform_on_sphere=False -> the original script's scheme: theta ~ Unif[0, pi],
                               which over-samples the poles (NOT area-uniform).
    """
    samples = []
    for _ in range(n_samples):
        phi = rng.uniform(0, 2 * np.pi, n_points)
        if uniform_on_sphere:
            theta = np.arccos(1.0 - 2.0 * rng.uniform(0.0, 1.0, n_points))
        else:
            theta = rng.uniform(0, np.pi, n_points)
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        samples.append(np.column_stack([x, y, z]))
    return samples


def first_zero_crossing(kappas, curve):
    """First lengthscale where `curve` crosses from >= 0 to < 0.

    Returns None if the curve never goes negative on the grid. Uses linear
    interpolation between the two bracketing grid points.
    """
    for i in range(len(curve) - 1):
        y0, y1 = curve[i], curve[i + 1]
        if y0 >= 0 and y1 < 0:
            x0, x1 = kappas[i], kappas[i + 1]
            # linear interpolation for the root
            return x0 + (x1 - x0) * y0 / (y0 - y1)
    return None


# --------------------------------------------------------------------------- #
# Core computation (cached)
# --------------------------------------------------------------------------- #
@st.cache_data(show_spinner=False)
def compute_curves(n_points, n_samples, kappa_min, kappa_max, n_kappa,
                   seed, uniform_on_sphere):
    """Compute per-sample minimum-eigenvalue curves and their zero-crossings."""
    rng = np.random.default_rng(seed)
    kappas = np.linspace(kappa_min, kappa_max, n_kappa)
    samples = generate_sphere_points(n_points, n_samples, rng, uniform_on_sphere)

    curves = np.empty((n_samples, n_kappa), dtype=float)
    crossings = []
    for s, points in enumerate(samples):
        D = build_squared_distance_matrix(points)
        for k, kap in enumerate(kappas):
            K = gaussian_gram(D, kap)
            curves[s, k] = np.linalg.eigvalsh(K).min()
        xc = first_zero_crossing(kappas, curves[s])
        if xc is not None:
            crossings.append(xc)

    mean_curve = curves.mean(axis=0)
    std_curve = curves.std(axis=0)
    mean_crossing_of_mean = first_zero_crossing(kappas, mean_curve)
    return {
        "kappas": kappas,
        "curves": curves,
        "mean_curve": mean_curve,
        "std_curve": std_curve,
        "crossings": np.array(crossings),
        "mean_crossing_of_mean": mean_crossing_of_mean,
    }


# --------------------------------------------------------------------------- #
# UI
# --------------------------------------------------------------------------- #
st.set_page_config(page_title="Geodesic RBF min-eigenvalue", layout="wide")
st.title("Minimum eigenvalue of the geodesic Gaussian kernel on $S^2$")
st.caption(
    "See `docs/eigenvalue_analysis.md` and `docs/schoenberg_intuition.md` for the "
    "theory behind the negative eigenvalues and the zero-crossing."
)

with st.sidebar:
    st.header("Parameters")
    n_points = st.slider("Points per sample (n)", 2, 200, 10)
    n_samples = st.slider("Number of samples", 1, 200, 20)
    st.subheader("Lengthscale grid")
    kappa_min = st.number_input("κ min", value=0.1, min_value=0.01, step=0.05)
    kappa_max = st.number_input("κ max", value=float(round(np.pi, 4)),
                                min_value=0.2, step=0.05)
    n_kappa = st.slider("κ resolution", 20, 400, 100)
    st.subheader("Sampling")
    uniform_on_sphere = st.checkbox(
        "Area-uniform sampling (recommended)", value=True,
        help="Off = the original script's theta~Unif[0,pi], which over-samples "
             "the poles and is NOT uniform on the sphere.",
    )
    seed = st.number_input("Random seed", value=0, min_value=0, step=1)

res = compute_curves(n_points, n_samples, float(kappa_min), float(kappa_max),
                     int(n_kappa), int(seed), bool(uniform_on_sphere))

kappas = res["kappas"]
curves = res["curves"]
mean_curve = res["mean_curve"]
std_curve = res["std_curve"]
crossings = res["crossings"]

# ----- Plot ----- #
fig, ax = plt.subplots(1, 1, figsize=(10, 7))
for curve in curves:
    ax.plot(kappas, curve, color="grey", alpha=0.3, linewidth=1)
ax.fill_between(kappas, mean_curve - std_curve, mean_curve + std_curve,
                color="lightgrey", alpha=0.5, label="±1 std")

colors = plt.cm.PuOr(kappas / kappas.max())
for i in range(len(kappas) - 1):
    ax.plot(kappas[i:i + 2], mean_curve[i:i + 2], color=colors[i],
            linewidth=3, alpha=0.9)

ax.axhline(y=0, color="black", linestyle="--", alpha=0.7, linewidth=2)

# Mark the mean of per-sample crossings
if len(crossings) > 0:
    ax.axvline(crossings.mean(), color="crimson", linestyle=":", linewidth=2,
               label=f"mean crossing κ* = {crossings.mean():.3f}")

mappable = plt.cm.ScalarMappable(cmap=plt.cm.PuOr,
                                 norm=plt.Normalize(vmin=0, vmax=kappas.max()))
mappable.set_array(kappas)
cbar = plt.colorbar(mappable, ax=ax, shrink=1.0, aspect=20)
cbar.set_label("lengthscale", fontsize=12)
cbar.ax.invert_yaxis()

ax.set_xlabel("length scale κ", fontsize=12)
ax.set_ylabel("Minimum Eigenvalue", fontsize=12)
ax.set_title("Minimum Eigenvalue of RBF Kernel", fontsize=12)
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_xlim([kappas[0], kappas[-1]])
plt.tight_layout()

col_plot, col_stats = st.columns([3, 1])
with col_plot:
    st.pyplot(fig)
with col_stats:
    st.subheader("Zero-crossing κ*")
    n_neg = len(crossings)
    st.metric("Samples that go negative", f"{n_neg} / {n_samples}")
    if n_neg > 0:
        st.metric("Mean crossing κ*", f"{crossings.mean():.4f}")
        st.metric("Std of crossing", f"{crossings.std():.4f}")
        st.metric("Min / Max", f"{crossings.min():.3f} / {crossings.max():.3f}")
    if res["mean_crossing_of_mean"] is not None:
        st.metric("Crossing of the MEAN curve",
                  f"{res['mean_crossing_of_mean']:.4f}")
    st.caption(
        "π/2 ≈ 1.5708 is the *mean pairwise distance*, a back-of-the-envelope "
        "scale — not the predicted crossing. The crossing is a spectral quantity "
        "and typically lands below π/2 (and drifts lower as n grows)."
    )
