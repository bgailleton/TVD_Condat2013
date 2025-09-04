import re
import sys
from pathlib import Path

increment = sys.argv[1]
if increment not in {"0.1", "0.01"}:
    raise SystemExit("increment must be 0.1 or 0.01")

setup_path = Path("setup.py")
conf_path = Path("docs/conf.py")

def read_version(text):
    m = re.search(r"__version__\s*=\s*['\"](\d+)\.(\d+)\.(\d+)['\"]", text)
    if not m:
        raise RuntimeError("version string not found")
    return [int(part) for part in m.groups()]

version_parts = read_version(setup_path.read_text())
if increment == "0.1":
    version_parts[1] += 1
    version_parts[2] = 0
else:
    version_parts[2] += 1
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
