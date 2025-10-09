# Release checklist

1. Update the version number in `setup.py` and `docs/conf.py` using the `bump_version.py`
   helper (or edit both files manually if preparing a one-off hotfix). Commit the change
   before starting the wheel build workflow.
2. Run `python -m build` locally to confirm the sdist and the platform wheel can be
   produced and pass `python wheelhouse/debug_wheel.py dist` to verify the archives.
3. Trigger the "Release wheels" workflow in GitHub Actions. The workflow invokes
   `cibuildwheel` on Linux, macOS, and Windows, validates the resulting archives with
   `wheelhouse/debug_wheel.py`, and publishes the release artefacts.
4. **Do not edit wheels in place.** Any post-processing that rewrites files inside a wheel
   (for example running `strip` or an antivirus scanner on `TVDCondat2013*.so/.pyd` inside
   the archive) changes the payload without updating the zip metadata. The debug script will
   report a CRC mismatch and the wheel will be rejected by `pip`. If you need a modified
   binary, rebuild the wheel from source instead of mutating an existing artefact.
