"""The vendored C++ core must still be the R package's.

This repository holds a copy of the core the R package owns, which is what
makes the two languages agree bit for bit. A copy can be edited locally by
accident, and then the agreement is gone with nothing to say so.

Two things guard it, and only the first lives here. This file checks the
vendored files against the manifest `tools/sync_core.py` wrote, so a local
edit is caught immediately. The real guard is the rest of the suite: if the
core stops reproducing R's fits, the fixture tests fail.
"""
import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
CORE = ROOT / "src" / "core"
MANIFEST = ROOT / "src" / "core_manifest.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def manifest():
    if not MANIFEST.exists():
        pytest.skip("no manifest; run tools/sync_core.py against a PACE checkout")
    return json.loads(MANIFEST.read_text())


def test_every_vendored_file_is_in_the_manifest(manifest):
    vendored = {path.name for path in CORE.glob("*.[ch]pp")}
    assert vendored == set(manifest["files"]), (
        "the vendored core and the manifest disagree on which files exist; "
        "re-run tools/sync_core.py")


def test_no_vendored_file_has_been_edited_locally(manifest):
    """An edit here would silently un-share the core with R."""
    drifted = [name for name, recorded in manifest["files"].items()
               if digest(CORE / name) != recorded]
    assert not drifted, (
        f"these vendored core files differ from the manifest: {drifted}. "
        "Change them in the PACE repository and re-sync, never here.")


def test_the_manifest_records_where_the_core_came_from(manifest):
    assert manifest["commit"], "the manifest should name the upstream commit"
    assert len(manifest["files"]) > 20
