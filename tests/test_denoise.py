import numpy as np

from onlpoisson.denoise import denoise_poisson


def test_denoise_shape_range_and_determinism():
    seed = 123
    rng1 = np.random.default_rng(seed)
    rng2 = np.random.default_rng(seed)

    x1 = rng1.poisson(lam=4.0, size=(8, 9)).astype(np.float64)
    x2 = rng2.poisson(lam=4.0, size=(8, 9)).astype(np.float64)

    y1 = denoise_poisson(x1, nw=1, npat=1, mu=0.2)
    y2 = denoise_poisson(x2, nw=1, npat=1, mu=0.2)

    assert y1.shape == x1.shape
    assert np.all(y1 >= 0.0)
    assert np.all(y1 <= 255.0)
    np.testing.assert_allclose(y1, y2, atol=0.0, rtol=0.0)
