"""
Streamlit app: the Euclidean-embeddability region on the sphere.

Place points on S^2 one at a time. At each step the sphere is shaded by

    lambda_max(P D P)   for the current set together with a candidate point x,

where D is the squared geodesic-distance matrix and P = I - 11^T/m projects onto
the zero-sum subspace. This value is >= 0 always, and equals 0 exactly when the
set (plus x) still embeds isometrically into Euclidean space -- i.e. when the
geodesic Gaussian kernel stays positive semidefinite for every lengthscale. The
"safe region" (where adding x keeps the set embeddable) is the near-zero band.

Once the current set is itself non-embeddable, adding any point keeps it so
(monotonicity), and no safe region remains -- which is the point of the demo.

Run with:
    uv run streamlit run embedding_region_app.py
"""

import numpy as np
import plotly.graph_objects as go
import streamlit as st

TOL = 1e-7


def sphere_grid(n_theta, n_phi):
    """Vertices of a (theta, phi) grid on the unit sphere, plus the meshgrid shape."""
    theta = np.linspace(0, np.pi, n_theta)
    phi = np.linspace(0, 2 * np.pi, n_phi)
    tt, pp = np.meshgrid(theta, phi, indexing="ij")
    x = np.sin(tt) * np.cos(pp)
    y = np.sin(tt) * np.sin(pp)
    z = np.cos(tt)
    V = np.stack([x, y, z], axis=-1).reshape(-1, 3)
    return V, (x, y, z), tt.shape


def geodesic(A, B):
    """Pairwise great-circle distances between rows of A and rows of B."""
    return np.arccos(np.clip(A @ B.T, -1.0, 1.0))


def lam_max_PDP(points):
    """lambda_max(P D P) for a set of points (>= 0; 0 iff embeddable)."""
    m = len(points)
    if m < 2:
        return 0.0
    D = geodesic(points, points) ** 2
    Pr = np.eye(m) - np.ones((m, m)) / m
    return float(np.linalg.eigvalsh(Pr @ D @ Pr).max())


@st.cache_data(show_spinner=False)
def embeddability_field(points_tuple, n_theta, n_phi):
    """Field lambda_max(P D P) over candidate next points x, on the sphere grid.

    Vectorized over all grid vertices via stacked eigvalsh.
    """
    pts = np.array(points_tuple, dtype=float).reshape(-1, 3)
    n = len(pts)
    V, (gx, gy, gz), shape = sphere_grid(n_theta, n_phi)
    G = len(V)
    m = n + 1

    # Base squared-distance block among the existing n points.
    base = geodesic(pts, pts) ** 2 if n else np.zeros((0, 0))
    # Squared distances from each grid vertex to the existing points: (G, n).
    cross = geodesic(V, pts) ** 2 if n else np.zeros((G, 0))

    M = np.zeros((G, m, m))
    if n:
        M[:, :n, :n] = base
        M[:, n, :n] = cross
        M[:, :n, n] = cross
    Pr = np.eye(m) - np.ones((m, m)) / m
    PMP = Pr @ M @ Pr  # broadcasts over the G axis
    field = np.linalg.eigvalsh(PMP).max(axis=-1).reshape(shape)
    return (gx, gy, gz), field


def find_safe_point(points_tuple, n_theta, n_phi):
    """Return a grid vertex that keeps the set embeddable, or None."""
    pts = np.array(points_tuple, dtype=float).reshape(-1, 3)
    V, _, _ = sphere_grid(n_theta, n_phi)
    n = len(pts)
    m = n + 1
    cross = geodesic(V, pts) ** 2 if n else np.zeros((len(V), 0))
    base = geodesic(pts, pts) ** 2 if n else np.zeros((0, 0))
    M = np.zeros((len(V), m, m))
    if n:
        M[:, :n, :n] = base
        M[:, n, :n] = cross
        M[:, :n, n] = cross
    Pr = np.eye(m) - np.ones((m, m)) / m
    field = np.linalg.eigvalsh(Pr @ M @ Pr).max(axis=-1)
    safe = np.where(field <= TOL)[0]
    return V[safe[0]] if len(safe) else None


# --------------------------------------------------------------------------- #
# State and controls
# --------------------------------------------------------------------------- #
st.set_page_config(page_title="Embeddability region on S^2", layout="wide")
st.title("Where can the next point go and keep the set Euclidean-embeddable?")
st.caption(
    "Sphere shaded by λ_max(P D P) for the current points plus a candidate x. "
    "0 (cool) = still embeddable; positive (warm) = adding x there breaks it. "
    "See `docs/embedding_and_curvature.md`."
)

if "pts" not in st.session_state:
    st.session_state.pts = []
if "rng_state" not in st.session_state:
    st.session_state.rng_state = 0

with st.sidebar:
    st.header("Build the configuration")
    res = st.select_slider("Grid resolution", options=[40, 60, 80, 120], value=60)
    n_theta, n_phi = res, 2 * res

    c1, c2 = st.columns(2)
    if c1.button("Add random point"):
        rng = np.random.default_rng(st.session_state.rng_state)
        st.session_state.rng_state += 1
        phi = rng.uniform(0, 2 * np.pi)
        theta = np.arccos(1 - 2 * rng.uniform(0, 1))
        st.session_state.pts.append(
            [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
        )
    if c2.button("Add safe point"):
        x = find_safe_point(tuple(map(tuple, st.session_state.pts)), n_theta, n_phi)
        if x is None:
            st.warning("No safe location remains.")
        else:
            st.session_state.pts.append(list(x))

    st.subheader("Add at latitude / longitude")
    lat = st.number_input("latitude (deg)", -90.0, 90.0, 0.0, step=5.0)
    lon = st.number_input("longitude (deg)", -180.0, 180.0, 0.0, step=5.0)
    if st.button("Add at lat/lon"):
        la, lo = np.radians(lat), np.radians(lon)
        st.session_state.pts.append(
            [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
        )

    if st.button("Reset"):
        st.session_state.pts = []

    st.metric("Points placed", len(st.session_state.pts))

# --------------------------------------------------------------------------- #
# Compute and render
# --------------------------------------------------------------------------- #
pts = st.session_state.pts
pts_tuple = tuple(map(tuple, pts))
base_lam = lam_max_PDP(np.array(pts).reshape(-1, 3)) if pts else 0.0
embeddable_now = base_lam <= TOL

(gx, gy, gz), field = embeddability_field(pts_tuple, n_theta, n_phi)
safe_frac = float((field <= TOL).mean())

col_plot, col_stats = st.columns([3, 1])

with col_plot:
    fig = go.Figure()
    fig.add_surface(
        x=gx, y=gy, z=gz, surfacecolor=field,
        colorscale="Turbo", cmin=0.0, cmax=max(field.max(), TOL),
        colorbar=dict(title="λ_max(PDP)"),
        showscale=True, opacity=1.0,
    )
    if pts:
        P = np.array(pts) * 1.02
        fig.add_scatter3d(
            x=P[:, 0], y=P[:, 1], z=P[:, 2], mode="markers",
            marker=dict(size=6, color="black"), name="points",
        )
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False),
            zaxis=dict(visible=False), aspectmode="data",
        ),
        margin=dict(l=0, r=0, t=0, b=0), height=650,
    )
    st.plotly_chart(fig, use_container_width=True)

with col_stats:
    st.subheader("Current set")
    st.metric("λ_max(PDP)", f"{base_lam:.4f}")
    if embeddable_now:
        st.success("Embeddable")
        st.metric("Safe region (next point)", f"{100 * safe_frac:.1f}% of sphere")
    else:
        st.error("Not embeddable")
        st.caption(
            "The set already fails to embed; adding any further point cannot "
            "repair it, so there is no safe region."
        )
