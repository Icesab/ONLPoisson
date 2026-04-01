"""Python implementation of the ONL denoiser from ``nlm_poisson004.c``."""

from __future__ import annotations

import math

import numpy as np


def _compute_cw(np_radius: int) -> np.ndarray:
    if np_radius <= 0:
        return np.array([1.0], dtype=np.float64)

    cw = np.zeros(np_radius + 1, dtype=np.float64)
    for i in range(1, np_radius + 1):
        s = 0.0
        for j in range(i, np_radius + 1):
            s += 1.0 / (np_radius * (2 * j + 1) * (2 * j + 1))
        cw[i] = s
    cw[0] = cw[1]
    return cw


def _mean_filter_reflect(img: np.ndarray, radius: int) -> np.ndarray:
    if radius <= 0:
        return img.astype(np.float64, copy=True)

    k = 2 * radius + 1
    padded = np.pad(img, radius, mode="reflect")
    integral = np.pad(np.cumsum(np.cumsum(padded, axis=0), axis=1), ((1, 0), (1, 0)))
    out = (
        integral[k:, k:]
        - integral[:-k, k:]
        - integral[k:, :-k]
        + integral[:-k, :-k]
    )
    return out / (k * k)


def denoise_onl(
    noisy: np.ndarray,
    nwin: int = 3,
    npat: int = 6,
    nlh: float = 0.5,
) -> np.ndarray:
    """Denoise a 2D Poisson-noisy image with ONL.

    Parameters mirror ``nlm_poisson004``:
    - ``nwin``: half of search-window size
    - ``npat``: half of patch size
    - ``nlh``: smoothing parameter used in exponential weights
    """

    noisy = np.asarray(noisy, dtype=np.float64)
    if noisy.ndim != 2:
        raise ValueError("`noisy` must be a 2D array.")
    if nwin < 0 or npat < 0:
        raise ValueError("`nwin` and `npat` must be non-negative.")
    if nlh <= 0:
        raise ValueError("`nlh` must be positive.")

    h, w = noisy.shape
    nwp = nwin + npat
    synoisy = np.pad(noisy, nwp, mode="reflect")

    avimage = np.sqrt(np.clip(_mean_filter_reflect(synoisy, nwin), 0.0, 255.0))

    # Restrict to original image support.
    avimage = avimage[nwp : nwp + h, nwp : nwp + w]

    patch_offsets = []
    patch_weights = []
    cw = _compute_cw(npat)
    for dx in range(-npat, npat + 1):
        for dy in range(-npat, npat + 1):
            patch_offsets.append((dx, dy))
            patch_weights.append(cw[max(abs(dx), abs(dy))])
    patch_weights = np.asarray(patch_weights, dtype=np.float64)

    out = np.zeros((h, w), dtype=np.float64)

    for i in range(h):
        x = i + nwp
        for j in range(w):
            y = j + nwp
            v = avimage[i, j]
            if v == 0:
                out[i, j] = 0.0
                continue

            center_patch = np.array(
                [synoisy[x + dx, y + dy] for dx, dy in patch_offsets], dtype=np.float64
            )

            num = 0.0
            den = 0.0
            for sx in range(x - nwin, x + nwin + 1):
                for sy in range(y - nwin, y + nwin + 1):
                    cand_patch = np.array(
                        [synoisy[sx + dx, sy + dy] for dx, dy in patch_offsets],
                        dtype=np.float64,
                    )
                    dist2 = np.sum(patch_weights * (center_patch - cand_patch) ** 2)
                    dist2 = max(math.sqrt(dist2) - 1.414 * v, 0.0)
                    wrho = math.exp(-dist2 / (nlh * v))
                    den += wrho
                    num += wrho * synoisy[sx, sy]

            out[i, j] = np.clip(num / den if den > 0 else synoisy[x, y], 0.0, 255.0)

    return out
