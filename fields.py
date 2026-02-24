"""Layer 1: Regulon field construction from transcript coordinates."""
import numpy as np
from scipy.ndimage import gaussian_filter


def make_grid(tx, ty, resolution=20.0, padding=50.0):
    """Create spatial grid covering transcript extent."""
    gx = np.arange(tx.min() - padding, tx.max() + padding, resolution)
    gy = np.arange(ty.min() - padding, ty.max() + padding, resolution)
    extent = (gx[0], gx[-1] + resolution, gy[0], gy[-1] + resolution)
    return gx, gy, extent


def transcript_kde(tx, ty, gx, gy, bandwidth=40.0, resolution=20.0):
    """Kernel density estimation for transcript positions."""
    nx, ny = len(gx), len(gy)
    ix = np.clip(((tx - gx[0]) / resolution).astype(int), 0, nx - 1)
    iy = np.clip(((ty - gy[0]) / resolution).astype(int), 0, ny - 1)
    counts = np.zeros((ny, nx))
    np.add.at(counts, (iy, ix), 1)
    return gaussian_filter(counts, sigma=bandwidth / resolution) / resolution**2


def build_gene_fields(df, genes, gx, gy, bandwidth=40.0, resolution=20.0):
    """Build KDE fields for all genes."""
    fields = {}
    for gene in genes:
        mask = df["feature_name"] == gene
        sub = df.loc[mask]
        fields[gene] = transcript_kde(
            sub["x_location"].values, sub["y_location"].values,
            gx, gy, bandwidth, resolution
        )
    return fields


def normalize_fields(fields):
    """Log-normalize gene density fields (analogous to scRNA-seq)."""
    genes = list(fields.keys())
    total = sum(fields[g] for g in genes)
    total = np.maximum(total, 1e-12)
    return {g: np.log(fields[g] / total * 1e4 + 1) for g in genes}, total


def zscore_fields(normalized_fields):
    """Z-score each gene field."""
    return {
        g: (f - f.mean()) / max(f.std(), 1e-12)
        for g, f in normalized_fields.items()
    }


def compute_regulon_field(zfields, targets, weights=None):
    """Weighted combination of target gene z-scored fields."""
    available = [g for g in targets if g in zfields]
    if not available:
        return None
    if weights is None:
        w = np.ones(len(available))
    else:
        w = np.array([weights[targets.index(g)] for g in available])
    result = sum(wi * zfields[g] for g, wi in zip(available, w))
    return result / max(np.sqrt(np.sum(w**2)), 1e-12)


def compute_gradient(field, resolution=20.0, smooth_sigma=1.5):
    """Gradient of a scalar field via Gaussian derivative filters."""
    f = gaussian_filter(field, sigma=smooth_sigma) if smooth_sigma > 0 else field
    gy, gx = np.gradient(f, resolution, resolution)
    return gx, gy, np.sqrt(gx**2 + gy**2)


def compute_laplacian(field, resolution=20.0, smooth_sigma=1.5):
    """Laplacian of a scalar field."""
    f = gaussian_filter(field, sigma=smooth_sigma) if smooth_sigma > 0 else field
    d2x = np.gradient(np.gradient(f, resolution, axis=1), resolution, axis=1)
    d2y = np.gradient(np.gradient(f, resolution, axis=0), resolution, axis=0)
    return d2x + d2y
