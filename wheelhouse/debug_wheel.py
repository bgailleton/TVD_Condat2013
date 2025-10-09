#!/usr/bin/env python3
"""Utilities for validating wheel archives produced by the CI."""

from __future__ import annotations

import sys
import zipfile
from pathlib import Path
from typing import Iterable

import zlib

HEADER_SIGNATURE = b"PK\x03\x04"


def _read_ignoring_crc(zf: zipfile.ZipFile, info: zipfile.ZipInfo) -> bytes:
    with zf.open(info) as handle:
        handle._update_crc = lambda *_args, **_kwargs: None  # type: ignore[attr-defined]
        return handle.read()


def _nearest_header_offset(path: Path, expected: int, search: int = 64) -> int | None:
    start = max(0, expected - search)
    with path.open("rb") as fh:
        fh.seek(start)
        blob = fh.read(search * 2)
    idx = blob.find(HEADER_SIGNATURE)
    if idx == -1:
        return None
    return start + idx


def _check_member(path: Path, zf: zipfile.ZipFile, info: zipfile.ZipInfo) -> Iterable[str]:
    try:
        with zf.open(info) as handle:
            for _ in iter(lambda: handle.read(1 << 14), b""):
                pass
    except zipfile.BadZipFile as exc:
        message = str(exc)
        details: list[str] = []
        if "Bad CRC-32" in message:
            payload = _read_ignoring_crc(zf, info)
            actual_crc = zlib.crc32(payload) & 0xFFFFFFFF
            details.append(
                "CRC mismatch: expected 0x%08x but archive contains 0x%08x"
                % (info.CRC, actual_crc)
            )
            details.append(
                "The member was modified after the wheel was built. Rebuild the wheel instead of editing it in place."
            )
        elif "Bad magic number" in message:
            actual = _nearest_header_offset(path, info.header_offset)
            if actual is None:
                details.append(
                    "Local file header missing at offset %d; the archive is truncated." % info.header_offset
                )
            else:
                delta = actual - info.header_offset
                details.append(
                    "Local file header shifted by %d bytes (expected %d, found %d)."
                    % (delta, info.header_offset, actual)
                )
                details.append(
                    "This indicates a post-processing step rewrote compressed data without updating zip metadata."
                )
        else:
            details.append(message)
        yield f"{info.filename}: {'; '.join(details)}"


def check_wheels(directory: str) -> int:
    """Validate wheels and return 0 if all are OK else 1."""

    status = 0
    for wheel in sorted(Path(directory).glob("*.whl")):
        try:
            with zipfile.ZipFile(wheel) as zf:
                issues = [msg for info in zf.infolist() for msg in _check_member(wheel, zf, info)]
                if issues:
                    status = 1
                    print(f"{wheel.name}: CORRUPTED")
                    for msg in issues:
                        print(f"  - {msg}")
                else:
                    print(f"{wheel.name}: OK")
        except zipfile.BadZipFile as exc:
            print(f"{wheel.name}: not a valid zip archive ({exc})")
            status = 1
    return status


if __name__ == "__main__":
    target_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    sys.exit(check_wheels(target_dir))
