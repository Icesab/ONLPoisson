"""Reference-oriented ONL Poisson denoising loop.

This module mirrors the core loop in ``nlm_poisson004.c`` with the requested
paper-aligned similarity update.
"""

from __future__ import annotations

import math

import numpy as np


NU = 1e-4


def _build_patch_weights(np_half: int, nxsy: int) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce ``w`` and ``dadr`` construction from the C implementation."""
    np2 = 2 * np_half + 1
    npnp = np2 * np2

    cw = np.zeros(np_half + 1, dtype=np.float64)
    for i in range(1, np_half + 1):
        acc = 0.0
        for j in range(i, np_half + 1):
            acc += 1.0 / (np_half * (2 * j + 1) * (2 * j + 1))
        cw[i] = acc
    if np_half >= 1:
        cw[0] = cw[1]
    else:
        cw[0] = 1.0

    w = np.empty(npnp, dtype=np.float64)
    dadr = np.empty(npnp, dtype=np.int64)
    idx = 0
    for x in range(-np_half, np_half + 1):
        for y in range(-np_half, np_half + 1):
            dadr[idx] = y * nxsy + x
            j = max(abs(x), abs(y))
            w[idx] = cw[j]
            idx += 1
    return w, dadr


def _avsqrt_image(synoisy: np.ndarray, nwp: int, nx: int, ny: int, nw: int, nxsy: int) -> np.ndarray:
    """Equivalent of ``avsqrt_image`` in C: local mean + clipping + sqrt."""
    nwnw = (2 * nw + 1) ** 2
    av = np.empty(nx * ny, dtype=np.float64)

    k = 0
    for x in range(nwp, nx + nwp):
        for y in range(nwp, ny + nwp):
            local = synoisy[y - nw : y + nw + 1, x - nw : x + nw + 1]
            mean_val = float(local.sum()) / nwnw
            mean_val = min(max(mean_val, 0.0), 255.0)
            av[k] = math.sqrt(mean_val)
            k += 1
    return av


def denoise_poisson(noisy: np.ndarray, nw: int = 3, npat: int = 6, mu: float = 0.5) -> np.ndarray:
    """Denoise a 2D image with ONL Poisson main loop.

    Args:
        noisy: Input image, 2D array-like.
        nw: Half search-window size.
        npat: Half patch size.
        mu: Similarity parameter ``mu`` (``nlh`` in C).

    Returns:
        Denoised image with same shape as ``noisy``.
    """
    noisy = np.asarray(noisy, dtype=np.float64)
    if noisy.ndim != 2:
        raise ValueError("noisy must be a 2D array")
    if nw < 0 or npat < 0:
        raise ValueError("nw and npat must be non-negative")

    nx, ny = noisy.shape
    nwp = nw + npat

    # Mirror/symmetric extension (C code's explicit symmetric padding intent).
    synoisy = np.pad(noisy, ((nwp, nwp), (nwp, nwp)), mode="symmetric")

    nxsy = nx + 2 * nwp
    avimage = _avsqrt_image(synoisy, nwp=nwp, nx=nx, ny=ny, nw=nw, nxsy=nxsy)

    w, dadr = _build_patch_weights(np_half=npat, nxsy=nxsy)
    npnp = w.size

    denoisy = np.empty(nx * ny, dtype=np.float64)
    syn_flat = synoisy.ravel(order="C")

    k = 0
    for x in range(nwp, nx + nwp):
        for y in range(nwp, ny + nwp):
            adr = y * nxsy + x
            v = avimage[k]
            if v == 0.0:
                denoisy[k] = 0.0
                k += 1
                continue

            uD = v * v
            weight_sum = 0.0
            value_sum = 0.0

            for xp in range(x - nw, x + nw + 1):
                for yp in range(y - nw, y + nw + 1):
                    dist2 = 0.0
                    adrp = yp * nxsy + xp
                    for j in range(npnp):
                        offset = int(dadr[j])
                        c1 = syn_flat[adr + offset]
                        c2 = syn_flat[adrp + offset]
                        diff = c1 - c2
                        dist2 += w[j] * diff * diff

                    dist2 = max(dist2 - 2.0 * uD, 0.0)
                    wrho = math.exp(-dist2 / (mu * v + NU))
                    weight_sum += wrho
                    value_sum += wrho * synoisy[yp, xp]

            out = value_sum / weight_sum
            denoisy[k] = min(max(out, 0.0), 255.0)
            k += 1

    return denoisy.reshape((nx, ny))
