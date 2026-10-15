"""Copy the C++ core in from a PACE checkout, and record what was copied.

The core is shared with the R package and is the reason the two languages
agree bit for bit. It lives in the R repository; this one vendors a copy.
Vendoring is the risk the split introduces, so it is made visible rather than
left to trust:

  * this script copies the files and writes core_manifest.json, a sha256 per
    file plus the upstream commit they came from;
  * tests/test_core_sync.py fails if the vendored files stop matching that
    manifest, which catches an accidental local edit;
  * the fixture tests fail if the core stops reproducing R's fits, which
    catches the thing that actually matters.

    python tools/sync_core.py ../PACE
"""
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
VENDORED = HERE / "src" / "core"
MANIFEST = HERE / "src" / "core_manifest.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def upstream_commit(checkout):
    try:
        return subprocess.run(["git", "-C", str(checkout), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    checkout = Path(argv[1]).resolve()
    source = checkout / "src" / "core"
    if not source.is_dir():
        print(f"no src/core under {checkout}; point this at a PACE checkout")
        return 1

    VENDORED.mkdir(parents=True, exist_ok=True)
    copied, changed = {}, []
    for path in sorted(source.glob("*.[ch]pp")):
        target = VENDORED / path.name
        before = digest(target) if target.exists() else None
        shutil.copy2(path, target)
        after = digest(target)
        copied[path.name] = after
        if before != after:
            changed.append(path.name)
    for stale in sorted(VENDORED.glob("*.[ch]pp")):
        if stale.name not in copied:
            stale.unlink()
            changed.append(f"{stale.name} (removed)")

    MANIFEST.write_text(json.dumps(
        {"source": str(checkout), "commit": upstream_commit(checkout), "files": copied},
        indent=2, sort_keys=True) + "\n")
    print(f"synced {len(copied)} files from {checkout} @ {upstream_commit(checkout)[:10]}")
    print("changed:", ", ".join(changed) if changed else "nothing")
    print("\nnow run the tests: a core that no longer reproduces R's fits will fail them.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
