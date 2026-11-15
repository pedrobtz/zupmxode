#!/bin/sh
# tools/ci/platform.sh <os>: zupmxode's per-OS checks, run once on each runner
# of platforms.yaml (os is macos, windows or ubuntu), with R_ACTIONS_PROFILE
# set to quick or full. Every check is a step here, in order; set -e names
# the one that failed. Run from the package root; it also runs locally.
#
#   1. src/core/ built and tested without R headers (tests/c), with R's
#      own C compiler.
#
# Later stages add theirs.
set -eu
os=${1:?usage: tools/ci/platform.sh <macos|windows|ubuntu>}
R=${R_HOME:+$R_HOME/bin/}R
command -v R >/dev/null 2>&1 && R=R

echo "::group::Toolchain"
cc=$("$R" CMD config CC)
echo "os = $os, CC = $cc"
$cc --version | head -1
echo "::endgroup::"

echo "::group::C core without R"
make -C tests/c clean check CC="$cc"
echo "::endgroup::"
