"""Opt-in Python startup hook used only for Shennong CLI subprocesses."""

from __future__ import annotations

import json
import os
from pathlib import Path


patches_json = os.environ.pop("SHENNONG_AUTOZYME_PATCHES", "")
status_file = os.environ.pop("SHENNONG_AUTOZYME_STATUS_FILE", "")
if patches_json:
    from shennong_autozyme import activate_autozyme

    try:
        patches = json.loads(patches_json)
    except (TypeError, ValueError):
        patches = []
    report = activate_autozyme(patches if isinstance(patches, list) else [])
    if status_file:
        try:
            Path(status_file).write_text(json.dumps(report, indent=2), encoding="utf-8")
        except OSError:
            pass
