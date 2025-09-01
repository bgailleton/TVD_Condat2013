"""Example illustrating 1D denoising with :mod:`TVDCondat2013`.

The script generates a piecewise-constant signal, corrupts it with Gaussian
noise, denoises it, and saves a figure comparing the three signals.
"""

import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path

from TVDCondat2013 import TVD


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
    lambda_tvd = 10.0
    denoised = TVD(noisy, lambda_tvd)

    # ------------------------------------------------------------------
    # Plot original, noisy, and denoised signals in separate panels
    # ------------------------------------------------------------------
    x = np.arange(clean.size)
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 6))

    axes[0].plot(x, clean, color="k", lw=1)
    axes[0].set_ylabel("Amplitude")
    axes[0].set_title("Original signal")

    axes[1].plot(x, noisy, color="0.6", lw=1)
    axes[1].set_ylabel("Amplitude")
    axes[1].set_title("Noisy signal")

    axes[2].plot(x, denoised, color="C1", lw=1.5)
    axes[2].set_xlabel("Sample")
    axes[2].set_ylabel("Amplitude")
    axes[2].set_title(f"Denoised with TVD (lambda={lambda_tvd})")

    fig.tight_layout()
    output_path = Path(__file__).with_name("Example.png")
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


if __name__ == "__main__":  # pragma: no cover - example script
    main()

