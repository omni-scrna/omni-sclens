#!/bin/sh
# Profile an entrypoint with denet, as an omnibenchmark entrypoint.
#
#   pca-prof: prof.sh pca.py
#
# Copied from rapids-singlecell/prof.sh. denet now comes from the env
# (almost-conductor, pinned in pixi.toml), not a host binary.
#
# Three things this works around, all verified against denet 0.6.0:
#
#  1. ob runs `./{entrypoint}` directly, so the conventional `denet pca.py`
#     form resolves to ./denet -- a file that does not exist in the module dir.
#     This wrapper IS the entrypoint, and `prof.sh pca.py` expands correctly.
#  2. denet needs a `run` subcommand and a `--` separator, or it parses the
#     module's own flags (--output_dir, --name) as its options.
#  3. denet cannot profile a shebang script: it execs the script, the kernel
#     swaps in the interpreter, and denet loses the process -- `cmd` comes back
#     empty and NO samples are written, while the script runs fine, so it fails
#     silently. The interpreter has to be named explicitly, so the shebang is
#     resolved here.
#
# FIFTH, denet >= 0.10: GPU sampling is OFF unless --gpu is passed. Without it
# the samples carry no VRAM at all, and nothing warns.
#
# FOURTH gotcha, found 2026-09-21: denet's sampling is ADAPTIVE. --interval is
# only the BASE; it backs off toward --max-interval (default 1000ms). Measured
# with --interval 50 alone: intervals ran 346, 233, 255, ... 586, 642ms, a mean
# of 376ms, so a 68ms `write` phase got ZERO samples and a 1669ms `pca` phase
# got one. Both flags must be pinned to get a fixed rate.
#
# Samples land beside the module's obkit-events.jsonl so the two align: denet
# gives RSS/threads/GPU every --interval ms, obkit gives phase boundaries, and
# the intersection is PER-PHASE memory. That is the only way to get compute
# memory uncontaminated by the load phase -- performance.txt reports one peak
# for the whole job, which for a module that materialises a matrix during load
# is the load's peak, not the algorithm's.
set -eu
# Pinned sampling rate; override per-run if a long job makes 50ms too many rows.
: "${DENET_INTERVAL_MS:=50}"
target="$1"
shift

out="."
prev=""
for a in "$@"; do
  [ "$prev" = "--output_dir" ] && out="$a"
  prev="$a"
done
mkdir -p "$out"

# Resolve the shebang: "#!/usr/bin/env python3" -> python3, "#!/usr/bin/Rscript"
# -> /usr/bin/Rscript. Falls back to executing directly if there is no shebang
# (a real binary, which denet handles).
shebang=$(head -1 "$target" | sed -n 's|^#! *||p')
case "$shebang" in
  */env\ *) interp=${shebang#*env } ;;
  "")       interp="" ;;
  *)        interp=$shebang ;;
esac

if [ -n "$interp" ]; then
  exec denet --json --out "$out/denet-samples.jsonl" \
       --interval "$DENET_INTERVAL_MS" --max-interval "$DENET_INTERVAL_MS" \
       --gpu --quiet run -- $interp "$target" "$@"
else
  exec denet --json --out "$out/denet-samples.jsonl" \
       --interval "$DENET_INTERVAL_MS" --max-interval "$DENET_INTERVAL_MS" \
       --gpu --quiet run -- "./$target" "$@"
fi
