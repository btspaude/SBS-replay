#!/usr/bin/env bash

# Request the raw CODA files needed for the first CDet LH2 timing-shift survey.
# The defaults cover the representative post-cross-target runs selected for
# the study.  Stream and segment limits are inclusive because CODA filenames
# start at zero.

set -euo pipefail

email="brash@jlab.org"
max_stream=2
max_segment=5
dry_run=0
runs=()

usage() {
  cat <<'EOF'
Usage: request_cdet_lh2_raw_to_cache.sh [options] [run ...]

Request GEp raw CODA files from tape with jcache. If no runs are supplied,
the representative LH2 survey runs 3649, 4346, 4724, 5295, and 5727 are used.

Options:
  --email ADDRESS       Completion-notification address (default: brash@jlab.org)
  --max-stream N        Highest stream number requested, inclusive (default: 2)
  --max-segment N       Highest segment number requested, inclusive (default: 5)
  --dry-run             Print the exact jcache command without submitting it
  -h, --help            Show this help

Examples:
  ./request_cdet_lh2_raw_to_cache.sh --dry-run
  ./request_cdet_lh2_raw_to_cache.sh 3649 4346
  ./request_cdet_lh2_raw_to_cache.sh --max-segment 9 5727
EOF
}

while (($#)); do
  case "$1" in
    --email)
      (($# >= 2)) || { echo "ERROR: --email requires an address" >&2; exit 2; }
      email=$2
      shift 2
      ;;
    --max-stream)
      (($# >= 2)) || { echo "ERROR: --max-stream requires a value" >&2; exit 2; }
      max_stream=$2
      shift 2
      ;;
    --max-segment)
      (($# >= 2)) || { echo "ERROR: --max-segment requires a value" >&2; exit 2; }
      max_segment=$2
      shift 2
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --*)
      echo "ERROR: unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      runs+=("$1")
      shift
      ;;
  esac
done

[[ $max_stream =~ ^[0-9]+$ ]] || { echo "ERROR: max stream must be a nonnegative integer" >&2; exit 2; }
[[ $max_segment =~ ^[0-9]+$ ]] || { echo "ERROR: max segment must be a nonnegative integer" >&2; exit 2; }

if ((${#runs[@]} == 0)); then
  runs=(3649 4346 4724 5295 5727)
fi

paths=()
for run in "${runs[@]}"; do
  [[ $run =~ ^[0-9]+$ ]] || { echo "ERROR: invalid run number: $run" >&2; exit 2; }
  for ((stream = 0; stream <= max_stream; ++stream)); do
    for ((segment = 0; segment <= max_segment; ++segment)); do
      paths+=("/mss/halla/sbs/GEp/raw/gep5_${run}.evio.${stream}.${segment}")
    done
  done
done

echo "Runs: ${runs[*]}"
echo "Streams: 0-${max_stream}; segments: 0-${max_segment}"
echo "Files requested: ${#paths[@]}"
echo "Notification: ${email}"

if ((dry_run)); then
  printf 'jcache get'
  printf ' %q' "${paths[@]}"
  printf ' -e %q\n' "$email"
  exit 0
fi

command -v jcache >/dev/null 2>&1 || {
  echo "ERROR: jcache is not available; run this script on a JLab farm host" >&2
  exit 127
}

jcache get "${paths[@]}" -e "$email"
