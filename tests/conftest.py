"""Pytest configuration ensuring the native extension is built before tests."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

try:
    import sysconfig
except ImportError as exc:  # pragma: no cover - Python always provides sysconfig
    raise RuntimeError("sysconfig module missing") from exc


def _build_extension() -> None:
    root = Path(__file__).resolve().parent.parent
    ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
    if ext_suffix is None:
        raise RuntimeError("Cannot determine extension suffix")

    target = root / f"TVDCondat2013{ext_suffix}"
    if target.exists():
        return

    subprocess.run(
        [sys.executable, "setup.py", "build_ext", "--inplace"],
        cwd=root,
        check=True,
    )


_build_extension()
