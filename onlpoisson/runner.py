"""Experiment runner reproducing the MATLAB ``demo_onl.m`` flow."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np

from .denoise import denoise_onl
from .metrics import fnmise

try:
    from scipy.io import loadmat
except Exception as exc:  # pragma: no cover
    loadmat = None
    _SCIPY_IMPORT_ERROR = exc
else:
    _SCIPY_IMPORT_ERROR = None


def _gaussian_kernel(size: int, sigma: float) -> np.ndarray:
    ax = np.arange(-(size // 2), size // 2 + 1)
    xx, yy = np.meshgrid(ax, ax, indexing="ij")
    k = np.exp(-(xx**2 + yy**2) / (2 * sigma * sigma))
    k /= k.sum()
    return k


def _imfilter_symmetric(img: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    kh, kw = kernel.shape
    ph, pw = kh // 2, kw // 2
    padded = np.pad(img, ((ph, ph), (pw, pw)), mode="reflect")
    out = np.empty_like(img, dtype=np.float64)
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            out[i, j] = np.sum(padded[i : i + kh, j : j + kw] * kernel)
    return out


def _load_images_mat(images_mat: str | Path) -> Dict[str, np.ndarray]:
    if loadmat is None:
        raise ImportError(
            "scipy is required to load images.mat"
        ) from _SCIPY_IMPORT_ERROR
    mat = loadmat(images_mat)
    return {
        "Spots": np.asarray(mat["spots"], dtype=np.float64),
        "Galaxy": np.asarray(mat["galaxy"], dtype=np.float64),
        "Ridges": np.asarray(mat["Ridges"], dtype=np.float64),
        "Barbara": np.asarray(mat["Barbara"], dtype=np.float64),
        "Cells": np.asarray(mat["cells"], dtype=np.float64),
    }


def run_table1_ours(
    images_mat: str | Path = "images.mat",
    repeats: int = 5,
    seed: int | None = 0,
) -> dict:
    """Run the Table-1 ONL reproduction loop and return per-image mean NMISE."""

    imgs = _load_images_mat(images_mat)
    names = ["Spots", "Galaxy", "Ridges", "Barbara", "Cells"]

    # Parameters copied from demo_onl.m (jjj=1..5 branches).
    params = {
        "Spots": dict(nlm0=9, nls0=6, nlwid=5, nlsi=1, nlh=0.15),
        "Galaxy": dict(nlm0=6, nls0=1, nlwid=5, nlsi=1, nlh=0.2),
        "Ridges": dict(nlm0=4, nls0=10, nlwid=7, nlsi=2, nlh=0.2),
        "Barbara": dict(nlm0=7, nls0=10, nlwid=1, nlsi=0, nlh=0.1),
        "Cells": dict(nlm0=3, nls0=6, nlwid=3, nlsi=1, nlh=0.2),
    }

    nmise_sum = {k: 0.0 for k in names}

    for rep in range(repeats):
        rng = np.random.default_rng(None if seed is None else seed + rep)
        for name in names:
            clean = imgs[name]
            p = params[name]

            noisy = rng.poisson(clean).astype(np.float64)
            den = denoise_onl(noisy, nwin=p["nlm0"], npat=p["nls0"], nlh=p["nlh"])

            if p["nlwid"] > 1:
                kernel = _gaussian_kernel(p["nlwid"], p["nlsi"])
                den = _imfilter_symmetric(den, kernel)

            nmise_sum[name] += fnmise(clean, den)

    nmise_mean = {k: nmise_sum[k] / repeats for k in names}
    return {
        "repeats": repeats,
        "seed": seed,
        "nmise_mean": nmise_mean,
    }
