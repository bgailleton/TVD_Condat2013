import re
import sys
from pathlib import Path

if len(sys.argv) != 2:
    raise SystemExit("usage: bump_version.py <increment>")

increment = sys.argv[1]

allowed_increments = {
    "+1": "major",
    "+0.1": "minor",
    "+0.0.1": "patch",
}

if increment not in allowed_increments:
    choices = ", ".join(sorted(allowed_increments))
    raise SystemExit(f"increment must be one of: {choices}")

setup_path = Path("setup.py")
conf_path = Path("docs/conf.py")

def read_version(text):
    m = re.search(r"__version__\s*=\s*['\"](\d+)\.(\d+)\.(\d+)['\"]", text)
    if not m:
        raise RuntimeError("version string not found")
    return [int(part) for part in m.groups()]

major, minor, patch = version_parts = read_version(setup_path.read_text())

kind = allowed_increments[increment]

if kind == "major":
    major += 1
    minor = 0
    patch = 0
elif kind == "minor":
    minor += 1
    patch = 0
else:
    patch += 1

version_parts = [major, minor, patch]
new_version = ".".join(str(p) for p in version_parts)

setup_text = re.sub(
    r"(__version__\s*=\s*['\"])([^'\"]+)(['\"])",
    lambda m: f"{m.group(1)}{new_version}{m.group(3)}",
    setup_path.read_text(),
)
setup_path.write_text(setup_text)

conf_text = conf_path.read_text()
conf_text = re.sub(
    r"(version = u')([^']+)('\n)",
    lambda m: f"{m.group(1)}{new_version}{m.group(3)}",
    conf_text,
)
conf_text = re.sub(
    r"(release = u')([^']+)('\n)",
    lambda m: f"{m.group(1)}{new_version}{m.group(3)}",
    conf_text,
)

conf_path.write_text(conf_text)

print(new_version)
