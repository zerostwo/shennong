"""Guarded AutoZyme activation helpers for Shennong Python backends."""

from __future__ import annotations

import importlib.metadata as metadata
import json
import os
from pathlib import Path
from typing import Iterable


AUTOZYME_SOURCE = (
    "github.com/zerostwo/autozyme@"
    "8fc2e9c3a7f70302f97589aaa9b0395dcf86f9bc#subdirectory=autozyme_py"
)

VALIDATED_UPSTREAM_VERSIONS = {
    "scanpy": {"scanpy": ("1.11.5",)},
    "scvelo": {"scvelo": ("0.3.5.dev1+gf63c0e705",)},
    "cell2location": {
        "cell2location": ("0.1.5",),
        "scvi-tools": ("1.3.3",),
        "pyro-ppl": ("1.9.1",),
    },
    "cellphonedb": {"cellphonedb": ("5.0.1",)},
}


def _requested_patches(patches: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(str(patch).strip() for patch in patches if str(patch).strip()))


def _base_report(patches: Iterable[str]) -> dict:
    requested = _requested_patches(patches)
    return {
        "provider": "autozyme",
        "source": AUTOZYME_SOURCE,
        "version": None,
        "requested": requested,
        "active": [],
        "inactive": requested.copy(),
        "disabled": bool(os.environ.get("AUTOZYME_DISABLED") or os.environ.get("AUTOZYME_DISABLE")),
        "upstream_versions": {},
        "compatibility_errors": {},
        "error": None,
    }


def _version_matches(actual: str, accepted: tuple[str, ...]) -> bool:
    return any(actual == expected or actual.startswith(f"{expected}+") for expected in accepted)


def _compatible_patches(report: dict) -> list[str]:
    compatible = []
    for patch in report["requested"]:
        problems = []
        for package, accepted in VALIDATED_UPSTREAM_VERSIONS.get(patch, {}).items():
            try:
                actual = metadata.version(package)
            except metadata.PackageNotFoundError:
                actual = None
            report["upstream_versions"][package] = actual
            if actual is None or not _version_matches(actual, accepted):
                problems.append(
                    f"{package} {actual or '(not installed)'}; requires "
                    f"{' or '.join(accepted)}"
                )
        if problems:
            report["compatibility_errors"][patch] = "; ".join(problems)
        else:
            compatible.append(patch)
    return compatible


def activate_autozyme(
    patches: Iterable[str], *, enabled: bool = True, skip_reason: str | None = None
) -> dict:
    """Activate requested patches and return evidence, never breaking the backend.

    AutoZyme itself owns parameter/version guards and delegates unsupported calls
    to the captured upstream implementation. This wrapper adds one more boundary:
    an unavailable or incompatible AutoZyme install is recorded, then the
    Shennong backend continues with the original Python package.
    """
    report = _base_report(patches)
    if not report["requested"]:
        return report
    if report["disabled"]:
        report["error"] = "AutoZyme activation was disabled by the environment."
        return report
    if not enabled:
        report["error"] = skip_reason or "AutoZyme activation was disabled for this call."
        return report
    try:
        import autozyme

        try:
            report["version"] = metadata.version("autozyme")
        except metadata.PackageNotFoundError:
            report["version"] = getattr(autozyme, "__version__", None)

        compatible = _compatible_patches(report)
        if not compatible:
            return report
        result = autozyme.activate(compatible)
        if isinstance(result, bool):
            result = {compatible[0]: result}
        status = autozyme.status()
        report["active"] = [
            patch
            for patch in report["requested"]
            if bool(result.get(patch)) and status.get(patch) == "active"
        ]
        report["inactive"] = [
            patch for patch in report["requested"] if patch not in report["active"]
        ]
    except Exception as error:  # AutoZyme is optional at the execution boundary.
        report["error"] = f"{type(error).__name__}: {error}"
    return report


def subprocess_autozyme_env(
    patches: Iterable[str], status_path: Path, base_env: dict[str, str] | None = None
) -> dict[str, str]:
    """Build an environment that activates AutoZyme in a Python CLI child."""
    env = dict(os.environ if base_env is None else base_env)
    shared_dir = str(Path(__file__).resolve().parent)
    current_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join(
        part for part in (shared_dir, current_pythonpath) if part
    )
    env["SHENNONG_AUTOZYME_PATCHES"] = json.dumps(_requested_patches(patches))
    env["SHENNONG_AUTOZYME_STATUS_FILE"] = str(status_path)
    return env


def read_autozyme_status(path: Path, patches: Iterable[str]) -> dict:
    """Read child-process activation evidence or return a guarded fallback."""
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(value, dict):
            return value
    except (OSError, ValueError, TypeError):
        pass
    report = _base_report(patches)
    report["error"] = "AutoZyme child activation did not produce a status record."
    return report
