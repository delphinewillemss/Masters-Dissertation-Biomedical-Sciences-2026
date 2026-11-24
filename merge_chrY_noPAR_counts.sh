#!/usr/bin/env bash
## uitvoeren met bash merge_chrY_noPAR_counts.sh (niet sbatch!)
set -Eeuo pipefail

BASE="/scratch/gent/vo/002/gvo00207/vsc52059"
RESULT_DIR="$BASE/chrY_noPAR_counts/sample_results"
SUMMARY="$BASE/chrY_noPAR_counts/chrY_noPAR_counts.tsv"

printf 'Sample\tReads_chrY_noPAR\n' > "$SUMMARY"
cat "$RESULT_DIR"/*.tsv >> "$SUMMARY"

echo "Merged $(ls "$RESULT_DIR"/*.tsv | wc -l) samples -> $SUMMARY"
ls "$BASE"/error/chrY_task_*.failed 2>/dev/null && echo "Let op: er zijn mislukte taken" || true