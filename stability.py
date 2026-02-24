"""Layer 3: Jacobian SVD, principal stretch, and Regulon Stability Index."""
import numpy as np
from .fields import compute_gradient


def compute_stability(regulon_fields, names, resolution=20.0):
    """Compute RSI and max principal stretch from Jacobian SVD.

    Parameters
    ----------
    regulon_fields : dict
        {name: 2D array} of regulon activity fields.
    names : list
        Ordered regulon names.
    resolution : float
        Grid spacing in microns.

    Returns
    -------
    sigma1 : ndarray (ny, nx)
        Maximum principal stretch (largest singular value of J).
    rsi : ndarray (ny, nx)
        Regulon Stability Index = sigma2 / sigma1.
    """
    ny, nx = regulon_fields[names[0]].shape
    grads = [(compute_gradient(regulon_fields[n], resolution)[0],
              compute_gradient(regulon_fields[n], resolution)[1]) for n in names]

    sigma1 = np.zeros((ny, nx))
    rsi = np.zeros((ny, nx))
    for y in range(ny):
        for x in range(nx):
            J = np.array([[gx[y, x], gy[y, x]] for gx, gy in grads])
            sv = np.linalg.svd(J, compute_uv=False)
            sigma1[y, x] = sv[0]
            rsi[y, x] = sv[1] / sv[0] if sv[0] > 1e-12 else 0.0

    return sigma1, rsi
