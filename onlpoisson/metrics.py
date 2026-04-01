"""Quality metrics for ONLPoisson."""

from __future__ import annotations

import numpy as np


def fnmise(original: np.ndarray, y_hat: np.ndarray) -> float:
    """Compute NMISE exactly as in ``fnmise.m``.

    Zero-valued pixels in ``original`` are discarded before computing
    ``mean((y_hat-original)^2 / original)``.
    """

    original = np.asarray(original, dtype=np.float64)
    y_hat = np.asarray(y_hat, dtype=np.float64)

    if original.shape != y_hat.shape:
        raise ValueError(
            f"Shape mismatch: original{original.shape} vs y_hat{y_hat.shape}."
        )

    vec_original = original.ravel()
    vec_y_hat = y_hat.ravel()

    mask = vec_original != 0
    if not np.any(mask):
        raise ValueError("NMISE is undefined when original is all zeros.")

    vec_original = vec_original[mask]
    vec_y_hat = vec_y_hat[mask]

    return float(np.mean((vec_y_hat - vec_original) ** 2 / vec_original))
