"""Example illustrating 1D denoising with :mod:`TVDCondat2013`.

The script generates a piecewise-constant signal, corrupts it with Gaussian
noise, denoises it with all available algorithms, and saves a figure comparing
the results.
"""

import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path

from TVDCondat2013 import fused_lasso, tvd_2013, tvd_2017, tvd_tautstring


def main() -> None:
    """Generate the figure demonstrating 1-D total-variation denoising."""
    # ------------------------------------------------------------------
    # Generate a piecewise-constant signal corrupted by Gaussian noise
    # ------------------------------------------------------------------
    rng = np.random.default_rng(0)
    segment_values = [0, 4, 1, 3, 0]
    segment_length = 100
    clean = np.repeat(segment_values, segment_length)
    noisy = clean + rng.normal(scale=1.0, size=clean.size)

    # Denoise the 1-D signal with a reasonable lambda
    lambda_tvd = 10.
    lambda_tvd2 = 1.5
    mu_l1 = 0.5
    denoised_v1 = tvd_2013(noisy, lambda_tvd)
    denoised_v2 = tvd_2017(noisy, lambda_tvd2)
    denoised_ts = tvd_tautstring(noisy, lambda_tvd)
    denoised_fl = fused_lasso(noisy, lambda_tvd, mu_l1)

    # ------------------------------------------------------------------
    # Plot original, noisy, and denoised signals in separate panels
    # ------------------------------------------------------------------
    x = np.arange(clean.size)
    fig, axes = plt.subplots(6, 1, sharex=True, figsize=(8, 12))

    axes[0].plot(x, clean, color="k", lw=1)
    axes[0].set_ylabel("Amplitude")
    axes[0].set_title("Original signal")

    axes[1].plot(x, noisy, color="0.6", lw=1)
    axes[1].set_ylabel("Amplitude")
    axes[1].set_title("Noisy signal")

    axes[2].plot(x, noisy, color="0.6", lw=1)
    axes[2].plot(x, denoised_v1, color="C1", lw=1.5)
    axes[2].set_ylabel("Amplitude")
    axes[2].set_title(f"tvd_2013 (lambda={lambda_tvd})")

    axes[3].plot(x, noisy, color="0.6", lw=1)
    axes[3].plot(x, denoised_v2, color="C2", lw=1.5)
    axes[3].set_ylabel("Amplitude")
    axes[3].set_title(f"tvd_2017 (lambda={lambda_tvd2})")

    axes[4].plot(x, noisy, color="0.6", lw=1)
    axes[4].plot(x, denoised_ts, color="C3", lw=1.5)
    axes[4].set_ylabel("Amplitude")
    axes[4].set_title(f"tvd_tautstring (lambda={lambda_tvd})")

    axes[5].plot(x, noisy, color="0.6", lw=1)
    axes[5].plot(x, denoised_fl, color="C4", lw=1.5)
    axes[5].set_xlabel("Sample")
    axes[5].set_ylabel("Amplitude")
    axes[5].set_title(
        f"fused_lasso (lambda={lambda_tvd}, mu={mu_l1})"
    )

    fig.tight_layout()
    output_path = Path(__file__).with_name("Example.png")
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


if __name__ == "__main__":  # pragma: no cover - example script
    main()

