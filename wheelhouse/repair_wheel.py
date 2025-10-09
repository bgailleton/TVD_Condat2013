#!/usr/bin/env python3
"""Repair wheels that were mutated in place by rebuilding their zip metadata."""

from __future__ import annotations

import argparse
import base64
import hashlib
import io
import struct
import zlib
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import zipfile

LOCAL_SIGNATURE = b"PK\x03\x04"
DATA_DESCRIPTOR_SIGNATURE = b"PK\x07\x08"


@dataclass
class Entry:
    name: str
    data: bytes
    compress_type: int
    date_time: tuple[int, int, int, int, int, int]


def _dos_to_tuple(dos_date: int, dos_time: int) -> tuple[int, int, int, int, int, int]:
    year = ((dos_date >> 9) & 0x7F) + 1980
    month = (dos_date >> 5) & 0x0F or 1
    day = dos_date & 0x1F or 1
    hour = (dos_time >> 11) & 0x1F
    minute = (dos_time >> 5) & 0x3F
    second = (dos_time & 0x1F) * 2
    return (year, month, day, hour, minute, second)


def _iter_local_records(blob: bytes) -> Iterable[Entry]:
    idx = 0
    length = len(blob)
    while idx < length:
        marker = blob.find(LOCAL_SIGNATURE, idx)
        if marker == -1:
            break
        if marker + 30 > length:
            break
        (ver, flag, comp, mtime, mdate, crc, csize, usize, name_len, extra_len) = struct.unpack_from(
            "<HHHHHIIIHH", blob, marker + 4
        )
        name_start = marker + 30
        name_end = name_start + name_len
        name = blob[name_start:name_end].decode("utf-8")
        extra_start = name_end
        data_start = extra_start + extra_len
        idx = data_start + csize

        if idx > length:
            break

        if flag & 0x08:
            descriptor = blob[idx:idx + 16]
            if descriptor.startswith(DATA_DESCRIPTOR_SIGNATURE):
                crc, csize, usize = struct.unpack_from("<III", descriptor, 4)
                idx += 16
            else:
                crc, csize, usize = struct.unpack_from("<III", blob, idx)
                idx += 12

        cdata = blob[data_start:data_start + csize]
        if comp == zipfile.ZIP_STORED:
            data = cdata
        elif comp == zipfile.ZIP_DEFLATED:
            data = zlib.decompress(cdata, -zlib.MAX_WBITS)
        else:
            raise RuntimeError(f"Unsupported compression method {comp} for {name}")

        if zlib.crc32(data) & 0xFFFFFFFF != crc:
            raise RuntimeError(f"Payload for {name} does not match recorded CRC")

        yield Entry(name=name, data=data, compress_type=comp, date_time=_dos_to_tuple(mdate, mtime))


def _rebuild_record(entries: OrderedDict[str, Entry], record_name: str) -> bytes:
    lines = []
    for name, entry in entries.items():
        if name == record_name:
            continue
        digest = base64.urlsafe_b64encode(hashlib.sha256(entry.data).digest()).rstrip(b"=").decode("ascii")
        lines.append(f"{name},sha256={digest},{len(entry.data)}")
    lines.append(f"{record_name},,")
    return ("\n".join(lines) + "\n").encode("utf-8")


def repair_wheel(path: Path, *, in_place: bool, output_dir: Path | None) -> Path:
    blob = path.read_bytes()
    entries = OrderedDict[str, Entry]()
    record_name: str | None = None

    for entry in _iter_local_records(blob):
        entries.pop(entry.name, None)
        entries[entry.name] = entry
        if entry.name.endswith(".dist-info/RECORD"):
            record_name = entry.name

    if not entries:
        raise RuntimeError(f"{path.name}: no local file records found")

    if record_name is None:
        raise RuntimeError(f"{path.name}: RECORD file missing")

    record_entry = entries.pop(record_name)
    rebuilt_record = _rebuild_record(entries, record_name)
    entries[record_name] = Entry(
        name=record_name,
        data=rebuilt_record,
        compress_type=record_entry.compress_type,
        date_time=record_entry.date_time,
    )

    if in_place and output_dir is not None:
        raise ValueError("--in-place and --output cannot be combined")

    if in_place:
        destination = path
    else:
        if output_dir is None:
            destination = path.with_suffix(".repaired.whl")
        else:
            output_dir.mkdir(parents=True, exist_ok=True)
            destination = output_dir / path.name

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as zf:
        for entry in entries.values():
            info = zipfile.ZipInfo(entry.name)
            info.compress_type = entry.compress_type
            info.date_time = entry.date_time
            zf.writestr(info, entry.data)

    buffer.seek(0)
    data = buffer.read()
    destination.write_bytes(data)
    return destination


def _expand_targets(targets: list[str]) -> list[Path]:
    paths: list[Path] = []
    for raw in targets:
        candidate = Path(raw)
        if candidate.is_dir():
            paths.extend(sorted(candidate.glob("*.whl")))
        else:
            paths.append(candidate)
    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("targets", nargs="+", help="Wheel files or directories containing wheels to repair")
    parser.add_argument("--in-place", action="store_true", help="Rewrite the wheels in place")
    parser.add_argument(
        "--output",
        type=Path,
        help="Directory where repaired wheels will be written (defaults to <wheel>.repaired.whl when omitted)",
    )
    args = parser.parse_args(argv)

    exit_code = 0
    for wheel in _expand_targets(args.targets):
        try:
            new_path = repair_wheel(wheel, in_place=args.in_place, output_dir=args.output)
        except Exception as exc:  # noqa: BLE001 - surface the failure
            print(f"{wheel}: {exc}")
            exit_code = 1
        else:
            if args.in_place:
                print(f"{wheel}: repaired in place")
            else:
                print(f"{wheel}: wrote {new_path}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
