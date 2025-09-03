#!/usr/bin/env python3
import sys
from pathlib import Path
import zipfile


def check_wheels(directory: str) -> None:
    """
    Validate every .whl file in `directory`. Prints whether each wheel is OK
    or identifies any corrupted member.
    """
    for wheel in Path(directory).glob("*.whl"):
        try:
            with zipfile.ZipFile(wheel) as zf:
                bad_entry = zf.testzip()
                if bad_entry:
                    print(f"{wheel.name}: corrupted member '{bad_entry}'")
                else:
                    print(f"{wheel.name}: OK")
        except zipfile.BadZipFile as exc:
            print(f"{wheel.name}: not a valid zip archive ({exc})")


if __name__ == "__main__":
    target_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    check_wheels(target_dir)
