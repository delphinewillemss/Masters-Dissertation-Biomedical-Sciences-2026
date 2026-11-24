##uitvoeren na sbatch chromosome_coverage.job met bash merge_chromosome_coverage.sh
## geeft een summary per sample en chromosoom weer met: Total_depth = Som van alle read-depthwaarden over alle posities, Mean_depth = Gemiddelde sequencing coverage, Covered_positions = Aantal posities met minstens één read, Positions_counted = Totaal aantal posities
#!/usr/bin/env bash

set -Eeuo pipefail

BASE="/scratch/gent/vo/002/gvo00207/vsc52059"

OUT_DIR="$BASE/chromosome_coverage"
RESULT_DIR="$OUT_DIR/sample_results"
ERROR_DIR="$BASE/error"

SUMMARY="$OUT_DIR/coverage_summary.tsv"
ERRORS="$OUT_DIR/processing_errors.tsv"

mkdir -p "$OUT_DIR" "$RESULT_DIR" "$ERROR_DIR"

# Create summary header
printf 'Sample\tChromosome\tTotal_depth\tMean_depth\tCovered_positions\tPositions_counted\n' > "$SUMMARY"

printf 'Sample\tProblem\n' > "$ERRORS"

shopt -s nullglob

# Combine all sample results
result_files=("$RESULT_DIR"/*.tsv)

for file in "${result_files[@]}"; do
    cat "$file" >> "$SUMMARY"
done

# Collect failed samples
failed_files=("$ERROR_DIR"/*.failed)

for file in "${failed_files[@]}"; do

    sample="$(basename "$file" .failed)"
    problem="$(tr '\t\n' '  ' < "$file")"

    printf '%s\t%s\n' "$sample" "$problem" >> "$ERRORS"

done

echo "Merged ${#result_files[@]} sample results."

echo "Summary: $SUMMARY"
echo "Errors: $ERRORS"