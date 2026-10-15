"""How closely the fixtures can be matched depends on how PACE was built.

The fixtures are the R package's own output, computed on an arm64 Mac with
Apple clang. What a Python build can reproduce depends on what changed, and
all three cases below are measured on CI rather than assumed:

  * same architecture and compiler -- bit-identical, every digit of U, alpha
    and the contamination loadings.
  * same architecture, different compiler (arm64 Linux, g++) -- every
    tolerance test passes; only exact equality is lost.
  * different architecture (x86-64) -- Eigen vectorises and reduces in a
    different order, and thirteen IRLS iterations amplify that to about
    2.7e-4 on 0.2% of the random effects.

So exact equality is asserted only on the reference build, and the numeric
tolerances gain a floor anywhere off the reference architecture. The floor is
not a fudge for a failing test: it is the measured size of the difference,
with room to spare.
"""
import os
import platform
import sys

MACHINE = platform.machine().lower()
ON_REFERENCE_ARCH = MACHINE in {"arm64", "aarch64"}

# The build the fixtures came from. Override to force either behaviour.
_forced = os.environ.get("PACE_REFERENCE_BUILD", "").strip().lower()
if _forced in {"1", "true", "yes"}:
    REFERENCE_BUILD = True
elif _forced in {"0", "false", "no"}:
    REFERENCE_BUILD = False
else:
    REFERENCE_BUILD = sys.platform == "darwin" and ON_REFERENCE_ARCH

# Measured: x86-64 differs from the arm64 fixtures by up to 2.7e-4. The floor
# is an order of magnitude above that, which still catches any real
# regression -- a genuine bug moves these numbers far further than this.
CROSS_ARCH_FLOOR = 0.0 if ON_REFERENCE_ARCH else 2e-3

WHERE = f"{sys.platform}/{MACHINE}"


def tol(atol, scale=1.0):
    """The tolerance for a comparison whose exact-build tolerance is `atol`.

    `scale` is the natural magnitude of what is being compared, so that one
    knob expresses the difference in each quantity's own units. The variance
    components are order one and take the default; a percentage runs to a
    hundred, and the same relative difference there is a hundred times
    larger in absolute terms. Without it the floor is tight by two orders of
    magnitude on exactly those columns.
    """
    return max(atol, CROSS_ARCH_FLOOR * scale)
