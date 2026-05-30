"""
Generate minimal eigenvalue plot vs lengthscale for different samples of points on the entire sphere.
Shows individual lines in grey, mean with colormap, and standard deviation in light grey background.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import math

def sph(lat_deg, lon_deg):
    """Unit vector on S^2 from latitude (deg) and longitude (deg)."""
    lat = math.radians(lat_deg)   # -90..90
    lon = math.radians(lon_deg)   # -180..180
    x = math.cos(lat)*math.cos(lon)
    y = math.cos(lat)*math.sin(lon)
    z = math.sin(lat)
    return np.array([x,y,z], dtype=float)

def geodesic_arc(u, v):
    """Great-circle distance (arc-length) on S^2 between unit vectors u, v."""
    dot = float(np.dot(u, v))
    dot = max(-1.0, min(1.0, dot))  # numerical safety
    angle = math.acos(dot)
    if angle > np.pi:
        angle = 2*np.pi - angle # Be sure to have the shortest path
    return angle

def gaussian_gram(Dsq, kappa):
    """Gaussian Gram from squared-distance matrix."""
    return np.exp(-Dsq / (2 * kappa**2))

def build_squared_distance_matrix(points_xyz, dist_fn):
    n = len(points_xyz)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            D[i, j] = dist_fn(points_xyz[i], points_xyz[j]) ** 2
    return D

def generate_sphere_points(n_points, n_samples):
    """Generate n_samples collections of n_points random points on the entire sphere."""
    all_samples = []
    
    for _ in range(n_samples):
        # Generate random points on entire sphere
        # Use spherical coordinates: theta in [0, pi], phi in [0, 2*pi]
        theta = np.random.uniform(0, np.pi, n_points)    # latitude from 0 to pi (entire sphere)
        phi = np.random.uniform(0, 2*np.pi, n_points)    # longitude from 0 to 2*pi
        
        # Convert to Cartesian coordinates
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        
        points = np.column_stack([x, y, z])
        all_samples.append(points)
    
    return all_samples

def create_eigenvalue_plot_with_samples():
    """Create eigenvalue plot for multiple samples of points on the entire sphere."""
    print("Creating eigenvalue plot with multiple samples...")
    
    # Parameters
    n_points = 10  # Number of points per sample
    n_samples = 20  # Number of different samples
    lambda_values = np.linspace(0.1, np.pi, 100)
    
    # Generate multiple samples of points on entire sphere
    all_samples = generate_sphere_points(n_points, n_samples)
    
    # Store all eigenvalue curves
    all_eigenval_curves = []
    
    # Compute eigenvalues for each sample
    for sample_idx, points in enumerate(all_samples):
        print(f"Processing sample {sample_idx + 1}/{n_samples}")
        
        # Build geodesic squared distance matrix for this sample
        D_geo_sq = build_squared_distance_matrix(points, geodesic_arc)
        
        # Store minimum eigenvalues for this sample
        min_eigenvals = []
        
        for lam in lambda_values:
            K = gaussian_gram(D_geo_sq, lam)
            eigenvals = np.linalg.eigvalsh(K)
            min_eigenval = eigenvals.min()
            min_eigenvals.append(min_eigenval)
        
        all_eigenval_curves.append(min_eigenvals)
    
    # Convert to numpy array for easier computation
    all_eigenval_curves = np.array(all_eigenval_curves)
    
    # Compute mean and standard deviation
    mean_eigenvals = np.mean(all_eigenval_curves, axis=0)
    std_eigenvals = np.std(all_eigenval_curves, axis=0)
    
    # Create the plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Plot individual curves in grey
    for curve in all_eigenval_curves:
        ax.plot(lambda_values, curve, color='grey', alpha=0.3, linewidth=1)
    
    # Plot standard deviation as light grey background
    ax.fill_between(lambda_values, 
                    mean_eigenvals - std_eigenvals, 
                    mean_eigenvals + std_eigenvals, 
                    color='lightgrey', alpha=0.5, label='±1 std')
    
    # Create colormap for mean curve
    colors = plt.cm.PuOr(lambda_values / np.pi)
    
    # Plot mean curve with colormap
    for i in range(len(lambda_values) - 1):
        ax.plot(lambda_values[i:i+2], mean_eigenvals[i:i+2], 
               color=colors[i], linewidth=3, alpha=0.9)
    
    # Add zero line
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.7, linewidth=2)
    
    # Add colorbar
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.PuOr, 
                                   norm=plt.Normalize(vmin=0, vmax=np.pi))
    mappable.set_array(lambda_values)
    cbar = plt.colorbar(mappable, ax=ax, shrink=1.0, aspect=20)
    cbar.set_label('lengthscale', fontsize=12)
    cbar.set_ticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    cbar.set_ticklabels(['0', 'π/4', 'π/2', '3π/4', 'π'])
    cbar.ax.invert_yaxis()
    
    # Formatting
    ax.set_xlabel('length scale', fontsize=12)
    ax.set_ylabel('Minimum Eigenvalue', fontsize=12)
    ax.set_title(f'Minimum Eigenvalue of RBF Kernel', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Set axis limits
    ax.set_xlim([lambda_values[0], lambda_values[-1]])
    
    plt.tight_layout()
    # plt.savefig('conflict/gifs/eigenvalue_samples_plot.png', dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
    
    print("Eigenvalue plot with samples saved as 'conflict/gifs/eigenvalue_samples_plot.png'")


def main():
    """Generate eigenvalue plot with multiple samples of points on the entire sphere."""
    print("Generating eigenvalue plot with multiple samples...")
    
    # Create output directory
    # os.makedirs('conflict/gifs', exist_ok=True)
    
    # Create the eigenvalue plot with samples
    create_eigenvalue_plot_with_samples()
    
    print("\nVisualization has been generated successfully!")
    print("- conflict/gifs/eigenvalue_samples_plot.png: Eigenvalue plot with multiple samples")

if __name__ == "__main__":
    main()
