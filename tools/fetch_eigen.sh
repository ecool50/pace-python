#!/usr/bin/env bash
# Fetch Eigen 3.4.0 into the project directory, for wheel builds.
#
# The version is pinned and the tarball checksummed because the fit reproduces
# the R package bit for bit only under Eigen 3.4: against 3.5 the random
# effects move by about 3e-4. A wheel built against the wrong Eigen would
# install cleanly and quietly disagree with the R package, which is the thing
# this whole build goes out of its way to prevent.
set -euo pipefail

VERSION=3.4.0
SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
URL="https://gitlab.com/libeigen/eigen/-/archive/${VERSION}/eigen-${VERSION}.tar.gz"
TARGET="eigen-${VERSION}"

if [ -d "${TARGET}" ]; then
  echo "eigen ${VERSION} already present at ${TARGET}"
  exit 0
fi

echo "fetching eigen ${VERSION}"
curl -fsSL "${URL}" -o eigen.tar.gz

# Compute and compare rather than piping to --check: macOS ships a BSD
# sha256sum whose flags differ from GNU's, so the usual one-liner fails there.
if command -v shasum >/dev/null 2>&1; then
  FOUND=$(shasum -a 256 eigen.tar.gz | cut -d" " -f1)
else
  FOUND=$(sha256sum eigen.tar.gz | cut -d" " -f1)
fi
if [ "${FOUND}" != "${SHA256}" ]; then
  echo "eigen tarball checksum mismatch" >&2
  echo "  expected ${SHA256}" >&2
  echo "  found    ${FOUND}" >&2
  exit 1
fi

tar xzf eigen.tar.gz
rm -f eigen.tar.gz
test -f "${TARGET}/Eigen/src/Core/util/Macros.h"
echo "eigen ${VERSION} ready at ${TARGET}"
