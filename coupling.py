"""Layer 2: Coupling tensor and phase boundary computation."""
import numpy as np
from scipy.ndimage import gaussian_filter


def compute_coupling_tensor(regulon_fields, names, delta=100.0, resolution=20.0):
    """Compute local covariance (coupling) tensor between regulon programs.

    Parameters
    ----------
    regulon_fields : dict
        {name: 2D array} of regulon activity fields.
    names : list
        Ordered regulon names.
    delta : float
        Coupling window radius in microns.
    resolution : float
        Grid spacing in microns.

    Returns
    -------
    coupling : ndarray (ny, nx, P, P)
        Local covariance tensor at each grid point.
    coupling_strength : ndarray (ny, nx)
        Frobenius norm ||C(s)||_F.
    eff_dim : ndarray (ny, nx)
        Effective dimensionality exp(-sum p_i log p_i).
    """
    P = len(names)
    ny, nx = regulon_fields[names[0]].shape
    stack = np.array([regulon_fields[n] for n in names])

    sigma_g = delta / resolution
    local_means = np.array([gaussian_filter(stack[i], sigma=sigma_g) for i in range(P)])
    centered = stack - local_means

    coupling = np.zeros((ny, nx, P, P))
    for i in range(P):
        for j in range(i, P):
            cov = gaussian_filter(centered[i] * centered[j], sigma=sigma_g)
            coupling[:, :, i, j] = cov
            coupling[:, :, j, i] = cov

    coupling_strength = np.zeros((ny, nx))
    eff_dim = np.zeros((ny, nx))
    for y in range(ny):
        for x in range(nx):
            C = coupling[y, x]
            coupling_strength[y, x] = np.sqrt(np.sum(C**2))
            eigvals = np.maximum(np.linalg.eigvalsh(C), 0)
            total = eigvals.sum()
            if total > 1e-12:
                p = eigvals[eigvals > 1e-12] / total
                eff_dim[y, x] = np.exp(-np.sum(p * np.log(p)))

    return coupling, coupling_strength, eff_dim


def compute_phase_boundaries(coupling, resolution=20.0, smooth_sigma=2.0):
    """Phase boundaries = ||grad(C)||: where coupling structure changes."""
    ny, nx, P, _ = coupling.shape
    pb = np.zeros((ny, nx))
    for i in range(P):
        for j in range(P):
            c_ij = gaussian_filter(coupling[:, :, i, j], sigma=smooth_sigma)
            gy, gx = np.gradient(c_ij, resolution, resolution)
            pb += gx**2 + gy**2
    return np.sqrt(pb)
