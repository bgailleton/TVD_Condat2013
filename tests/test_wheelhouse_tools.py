from __future__ import annotations

import subprocess
import sys
import zipfile
from pathlib import Path


def _copy_wheel(tmp_path: Path, filename: str) -> Path:
    src = Path(__file__).resolve().parent.parent / "wheelhouse" / filename
    dst = tmp_path / filename
    dst.write_bytes(src.read_bytes())
    return dst


def test_debug_reports_duplicate_eocd(tmp_path: Path) -> None:
    wheel = _copy_wheel(tmp_path, "tvdcondat2013-0.0.1-cp310-cp310-win32.whl")
    proc = subprocess.run(
        [sys.executable, "wheelhouse/debug_wheel.py", str(tmp_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "Multiple end-of-central-directory records" in proc.stdout
    assert "CORRUPTED" in proc.stdout


def test_repair_wheel_rewrites_metadata(tmp_path: Path) -> None:
    wheel = _copy_wheel(tmp_path, "tvdcondat2013-0.0.1-cp310-cp310-win32.whl")

    subprocess.run(
        [sys.executable, "wheelhouse/repair_wheel.py", "--in-place", str(wheel)],
        check=True,
    )

    proc = subprocess.run(
        [sys.executable, "wheelhouse/debug_wheel.py", str(tmp_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "OK" in proc.stdout and "CORRUPTED" not in proc.stdout

    with zipfile.ZipFile(wheel) as zf:
        record = zf.read("tvdcondat2013-0.0.1.dist-info/RECORD").decode("utf-8")
    assert "tvdcondat2013-0.0.1.dist-info/licenses/LICENSE,sha256=" in record
    assert ",1418" in record.splitlines()[1]
