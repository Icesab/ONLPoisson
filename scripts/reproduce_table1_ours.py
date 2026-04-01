#!/usr/bin/env python3
"""Reproduce ONL 'Ours' pipeline for Table 1.

Workflow (per image, per peak, 30 trials):
1) u = peak * img / img.max()
2) z ~ Poisson(u)
3) u1 = denoise_onl(z, d=13, D=11, mu=0.2, nu=1e-4)
4) u2 = gaussian_filter(u1, kernel=5, sigma=1, mode='reflect')
5) nmise = fnmise(u, u2)
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np

try:
    from scipy.io import loadmat  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    loadmat = None


IMAGE_NAMES = ["Spots", "Galaxy", "Cells", "Texture", "Fermi"]
PEAKS = [0.5, 1, 2, 3, 4, 5]
NUM_TRIALS = 30


def fnmise(u: np.ndarray, uhat: np.ndarray) -> float:
    """MATLAB fnmise equivalent: average ((uhat-u)^2 / u) over u != 0."""
    u_vec = np.asarray(u, dtype=np.float64).ravel()
    uhat_vec = np.asarray(uhat, dtype=np.float64).ravel()
    mask = u_vec != 0
    if not np.any(mask):
        return 0.0
    return float(np.mean(((uhat_vec[mask] - u_vec[mask]) ** 2) / u_vec[mask]))


def gaussian_kernel(size: int = 5, sigma: float = 1.0) -> np.ndarray:
    if size % 2 != 1:
        raise ValueError("Gaussian kernel size must be odd.")
    radius = size // 2
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    xx, yy = np.meshgrid(x, x, indexing="xy")
    k = np.exp(-(xx * xx + yy * yy) / (2.0 * sigma * sigma))
    k /= np.sum(k)
    return k


def reflect_pad(img: np.ndarray, pad: int) -> np.ndarray:
    return np.pad(img, ((pad, pad), (pad, pad)), mode="reflect")


def conv2_reflect(img: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    k = np.asarray(kernel, dtype=np.float64)
    if k.ndim != 2 or k.shape[0] != k.shape[1] or k.shape[0] % 2 != 1:
        raise ValueError("kernel must be odd square 2D array")
    pad = k.shape[0] // 2
    padded = reflect_pad(np.asarray(img, dtype=np.float64), pad)
    out = np.empty_like(img, dtype=np.float64)
    kh, kw = k.shape
    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            patch = padded[i : i + kh, j : j + kw]
            out[i, j] = np.sum(patch * k)
    return out


def denoise_onl(z: np.ndarray, d: int = 13, D: int = 11, mu: float = 0.2, nu: float = 1e-4) -> np.ndarray:
    """Pure NumPy ONL-style denoiser.

    Parameter mapping:
    - d: patch size (odd), patch radius = d//2
    - D: search size (odd), search radius = D//2
    - mu, nu: exponential weight denominator mu * sqrt(u_D) + nu
    """
    if d % 2 != 1 or D % 2 != 1:
        raise ValueError("d and D must be odd")

    z = np.asarray(z, dtype=np.float64)
    h, w = z.shape
    pr = d // 2
    sr = D // 2

    # Weighted patch metric (Gaussian weights, normalized)
    pw = gaussian_kernel(size=d, sigma=max(d / 6.0, 1e-6))

    # u_D estimate via local mean with (D x D) uniform window
    local_mean_kernel = np.full((D, D), 1.0 / (D * D), dtype=np.float64)
    u_D = conv2_reflect(z, local_mean_kernel)
    v = np.sqrt(np.clip(u_D, 0.0, None))

    zpad = reflect_pad(z, pr + sr)
    out = np.empty_like(z, dtype=np.float64)

    for i in range(h):
        for j in range(w):
            ii = i + pr + sr
            jj = j + pr + sr
            ref_patch = zpad[ii - pr : ii + pr + 1, jj - pr : jj + pr + 1]

            denom = mu * v[i, j] + nu
            wsum = 0.0
            acc = 0.0

            for di in range(-sr, sr + 1):
                for dj in range(-sr, sr + 1):
                    ci = ii + di
                    cj = jj + dj
                    nbr_patch = zpad[ci - pr : ci + pr + 1, cj - pr : cj + pr + 1]
                    dist2 = float(np.sum(pw * (ref_patch - nbr_patch) ** 2))
                    dist2 = max(dist2 - 2.0 * u_D[i, j], 0.0)
                    wrho = float(np.exp(-dist2 / denom))
                    wsum += wrho
                    acc += wrho * zpad[ci, cj]

            out[i, j] = acc / wsum if wsum > 0 else z[i, j]

    return out


def gaussian_filter(img: np.ndarray, kernel: int = 5, sigma: float = 1.0, mode: str = "reflect") -> np.ndarray:
    if mode != "reflect":
        raise ValueError("Only mode='reflect' is supported.")
    return conv2_reflect(np.asarray(img, dtype=np.float64), gaussian_kernel(kernel, sigma))


def _extract_image_value(raw: object) -> np.ndarray:
    arr = np.asarray(raw)
    while arr.dtype == object and arr.size == 1:
        arr = np.asarray(arr.item())
    return np.asarray(arr, dtype=np.float64)


def load_images(images_mat: Path, data_dir: Path) -> Dict[str, np.ndarray]:
    """Load Spots/Galaxy/Cells/Texture/Fermi from images.mat, fallback to data/ for missing images."""
    out: Dict[str, np.ndarray] = {}

    if images_mat.exists() and loadmat is not None:
        mat = loadmat(images_mat)
        key_variants = {
            "Spots": ["spots", "Spots"],
            "Galaxy": ["galaxy", "Galaxy"],
            "Cells": ["cells", "Cells"],
            "Texture": ["texture", "Texture"],
            "Fermi": ["fermi", "Fermi"],
        }
        for name, keys in key_variants.items():
            for k in keys:
                if k in mat:
                    out[name] = _extract_image_value(mat[k])
                    break

    missing = [name for name in IMAGE_NAMES if name not in out]
    if missing:
        for name in missing:
            npy_path = data_dir / f"{name}.npy"
            if npy_path.exists():
                out[name] = np.load(npy_path).astype(np.float64)

    still_missing = [name for name in IMAGE_NAMES if name not in out]
    if still_missing:
        raise FileNotFoundError(
            "Missing required images: "
            + ", ".join(still_missing)
            + ". Provide them in images.mat or as data/<Name>.npy"
        )

    return out


def run_once(img: np.ndarray, peak: float, rng: np.random.Generator) -> float:
    img = np.asarray(img, dtype=np.float64)
    imax = float(np.max(img))
    if imax <= 0:
        raise ValueError("Image max must be > 0")

    u = peak * img / imax
    z = rng.poisson(u).astype(np.float64)
    u1 = denoise_onl(z, d=13, D=11, mu=0.2, nu=1e-4)
    u2 = gaussian_filter(u1, kernel=5, sigma=1, mode="reflect")
    return fnmise(u, u2)


def summarize_rows(rows: Iterable[Dict[str, object]]) -> str:
    lines = ["Table1 Ours reproduction summary (NMISE):", "image, peak, mean, std"]
    for r in rows:
        lines.append(f"{r['image']}, {r['peak']}, {r['nmise_mean']:.6f}, {r['nmise_std']:.6f}")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--images-mat", default="images.mat", type=Path)
    parser.add_argument("--data-dir", default="data", type=Path)
    parser.add_argument("--output", default="results/table1_ours.csv", type=Path)
    parser.add_argument("--trials", default=NUM_TRIALS, type=int)
    args = parser.parse_args()

    images = load_images(args.images_mat, args.data_dir)
    rng = np.random.default_rng(0)

    rows: List[Dict[str, object]] = []

    for image_name in IMAGE_NAMES:
        img = images[image_name]
        for peak in PEAKS:
            vals = [run_once(img, peak, rng) for _ in range(args.trials)]
            rows.append(
                {
                    "image": image_name,
                    "peak": peak,
                    "trials": args.trials,
                    "nmise_mean": float(np.mean(vals)),
                    "nmise_std": float(np.std(vals, ddof=0)),
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["image", "peak", "trials", "nmise_mean", "nmise_std"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(summarize_rows(rows))
    print(f"\nSaved CSV: {args.output}")


if __name__ == "__main__":
    main()
