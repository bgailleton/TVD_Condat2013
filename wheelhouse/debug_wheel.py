#!/usr/bin/env python3
import sys
from pathlib import Path
import zipfile


def check_wheels(directory: str) -> int:
    """Validate wheels and return 0 if all are OK else 1."""
    status = 0
    for wheel in Path(directory).glob("*.whl"):
        try:
            with zipfile.ZipFile(wheel) as zf:
                bad_entry = zf.testzip()
                if bad_entry:
                    print(f"{wheel.name}: corrupted member '{bad_entry}'")
                    status = 1
                else:
                    print(f"{wheel.name}: OK")
        except zipfile.BadZipFile as exc:
            print(f"{wheel.name}: not a valid zip archive ({exc})")
            status = 1
    return status


if __name__ == "__main__":
    target_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    sys.exit(check_wheels(target_dir))
