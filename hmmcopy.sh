#!/bin/bash
# Gebruik:  bash hmmcopy.sh (niet indienen met sbatch!)
#   stappen: 1=split  2=map  3=mergemap  4=hmmcopy  5=merge  6=plots

BASE=/scratch/gent/vo/002/gvo00207/vsc52059
WORK=$BASE/hmmcopy
LOGS=$WORK/logs
REF=$BASE/mm10_ref/mm10.fa
SPLIT=$WORK/reference/mapsplit
UTILS=$BASE/hmmcopy_utils/util/mappability/internal
CHRLIST=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
BIN=500000
READLEN=100
STEP=${STEP:-submit}

case "$STEP" in

# INDIENEN (dit gebeurt als je "bash hmmcopy.sh" typt)
submit)
  SELF=$(readlink -f "$0")
  mkdir -p $LOGS $WORK/results/plots

  if [ -n "$1" ]; then START=$1
  elif [ -s $WORK/reference/mm10.map.wig ]; then
    START=4; echo "mm10.map.wig bestaat al -> start bij stap 4"
  else START=1; fi

  LAST=""
  submit_job() {   # naam  afterok|afterany  logbestand  [sbatch-opties]
    local name=$1 dtype=$2 out=$3; shift 3
    local dep=""; [ -n "$LAST" ] && dep="--dependency=$dtype:$LAST"
    local id
    id=$(sbatch --parsable $dep --job-name=$name --output=$LOGS/$out \
         --export=ALL,STEP=$name,SELF=$SELF "$@" "$SELF") \
      || { echo "FOUT bij indienen van $name"; exit 1; }
    LAST=${id%%;*}
    echo "$name -> job $LAST"
  }

  echo "Start bij stap $START"
  [ $START -le 1 ] && submit_job split     afterok  "s1_split_%j.out"        --cpus-per-task=1 --mem=4G --time=00:30:00
  [ $START -le 2 ] && submit_job map       afterok  "s2_map_%A_%a.out"       --cpus-per-task=1 --mem=8G --time=02:00:00 --array=1-21
  [ $START -le 3 ] && submit_job mergemap  afterok  "s3_mergemap_%j.out"     --cpus-per-task=1 --mem=4G --time=00:20:00
  [ $START -le 4 ] && submit_job hmmcopy   afterok  "s4_hmmcopy_%A_%a.out"   --cpus-per-task=1 --mem=8G --time=01:00:00 --array=1-86%20
  [ $START -le 5 ] && submit_job merge     afterany "s5_merge_%j.out"        --cpus-per-task=1 --mem=2G --time=00:10:00
  [ $START -le 6 ] && submit_job plots     afterany "s6_plots_%j.out"        --cpus-per-task=1 --mem=4G --time=00:30:00
  echo "Alles ingediend. Volgen met: squeue -u \$USER"
  exit 0
  ;;

# STAP 1: genoom splitsen per chromosoom
split)
  set -eo pipefail
  for e in 1 2 3 4 rev.1 rev.2; do
    [ -s ${REF}.$e.ebwt ] || { echo "FOUT: bowtie-index ${REF}.$e.ebwt ontbreekt"; exit 1; }
  done
  mkdir -p $SPLIT
  module purge
  module load SAMtools/1.22.1-GCC-14.2.0
  samtools faidx $REF
  for c in "${CHRLIST[@]}"; do samtools faidx $REF $c > $SPLIT/$c.fa; done
  echo "Klaar"
  ;;

# STAP 2: mappability per chromosoom (array 1-21)
map)
  set -eo pipefail
  CHR=${CHRLIST[$((SLURM_ARRAY_TASK_ID-1))]}
  LEN=$(awk -v c=$CHR '$1==c{print $2}' ${REF}.fai)
  [ -n "$LEN" ] || { echo "FOUT: lengte van $CHR niet gevonden"; exit 1; }
  OUT=$SPLIT/$CHR.map.wig
  [ -s $OUT ] && { echo "$CHR bestaat al, overgeslagen"; exit 0; }
  module purge
  module load Bowtie/1.3.1-GCC-13.3.0
  echo "$CHR (lengte $LEN) gestart: $(date)"
  $UTILS/fastaToRead -w $READLEN $SPLIT/$CHR.fa \
    | bowtie $REF -f - -v 0 --quiet \
    | perl $UTILS/readToMap.pl --max_hits 4 \
    | awk -v chr=$CHR -v len=$LEN -v bin=$BIN '
        /^fixedStep/ {next}
        { b=int(n/bin); s[b]+=$1; n++ }
        END {
          nb=int((len+bin-1)/bin)
          printf "fixedStep chrom=%s start=1 step=%d span=%d\n", chr, bin, bin
          for (i=0; i<nb; i++) {
            bl=bin; if ((i+1)*bin>len) bl=len-i*bin
            printf "%.6f\n", s[i]/bl
          }
        }' > $OUT.tmp
  NB=$(( (LEN + BIN - 1) / BIN ))
  [ $(wc -l < $OUT.tmp) -eq $((NB+1)) ] || { echo "FOUT: verkeerd aantal bins voor $CHR"; exit 1; }
  mv $OUT.tmp $OUT
  echo "$CHR klaar: $(date)"
  ;;

# STAP 3: mappability samenvoegen en controleren
mergemap)
  set -eo pipefail
  cd $WORK/reference
  : > mm10.map.wig.tmp
  for c in "${CHRLIST[@]}"; do cat $SPLIT/$c.map.wig >> mm10.map.wig.tmp; done
  diff <(grep fixedStep mm10.gc.wig) <(grep fixedStep mm10.map.wig.tmp) > /dev/null \
    || { echo "FOUT: chromosomen/bins komen niet overeen met mm10.gc.wig"; exit 1; }
  [ $(wc -l < mm10.gc.wig) -eq $(wc -l < mm10.map.wig.tmp) ] \
    || { echo "FOUT: aantal regels verschilt van mm10.gc.wig"; exit 1; }
  mv mm10.map.wig.tmp mm10.map.wig
  echo "mm10.map.wig OK: $(wc -l < mm10.map.wig) regels"
  grep -v fixedStep mm10.map.wig | awk '{s+=$1;n++} END{print "Gemiddelde mappability:", s/n, "| bins:", n}'
  rm -f $SPLIT/*.fa $SPLIT/*.fa.fai
  ;;

# STAP 4: HMMcopy per sample (array 1-86)
hmmcopy)
  set -eo pipefail
  cd $WORK
  WIG=$(ls readcounts/*.wig | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
  [ -n "$WIG" ] || { echo "Geen wig voor taak ${SLURM_ARRAY_TASK_ID}"; exit 1; }
  echo "Taak ${SLURM_ARRAY_TASK_ID}: $WIG"
  RF=$(mktemp --suffix=.R); trap 'rm -f $RF' EXIT
  cat > $RF << 'RSCRIPT'
library(HMMcopy)
args <- commandArgs(trailingOnly = TRUE)
f    <- args[1]
gc   <- "reference/mm10.gc.wig"
map  <- "reference/mm10.map.wig"
outdir <- "results"
chrs <- c(paste0("chr", 1:19), "chrX", "chrY")
state_labels <- c("1_HOMD", "2_HETD", "3_NEUT", "4_GAIN", "5_AMPL", "6_HLAMP")
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
s    <- sub("\\.wig$", "", basename(f))
outf <- file.path(outdir, paste0(s, "_segments.txt"))
if (file.exists(outf)) { message("Overslaan (bestaat al): ", s); quit(status = 0) }
message("Verwerk ", s)
unc  <- wigsToRangedData(f, gc, map)
corr <- correctReadcount(unc)
corr$chr <- factor(corr$chr, levels = chrs)
seg  <- HMMsegment(corr, verbose = FALSE)
segs <- seg$segs
segs$chr <- factor(segs$chr, levels = chrs)
segs <- segs[order(segs$chr, segs$start), ]
segs$state_label <- state_labels[segs$state]
saveRDS(list(corrected = corr, segs = seg), file.path(outdir, paste0(s, ".rds")))
write.table(segs, outf, sep = "\t", quote = FALSE, row.names = FALSE)
tryCatch({
  pdf(file.path(outdir, "plots", paste0(s, ".pdf")), width = 14, height = 5)
  plotSegments(corr, seg, main = s)
  dev.off()
}, error = function(e) { try(dev.off(), silent = TRUE)
                         message("plot mislukt voor ", s, ": ", conditionMessage(e)) })
message("Klaar: ", s)
RSCRIPT
  module purge
  module load R-bundle-Bioconductor/3.23-foss-2026.1-R-4.6.1
  export OPENBLAS_NUM_THREADS=1
  Rscript $RF "$WIG"
  ;;

# STAP 5: alle segmenten in één tabel
merge)
  cd $WORK/results
  first=$(ls *_segments.txt | head -1)
  { head -n 1 "$first" | awk '{print $0"\tsample"}'
    for f in *_segments.txt; do
      s=${f%_segments.txt}
      tail -n +2 "$f" | awk -v s="$s" '{print $0"\t"s}'
    done; } > ALL_segments.txt
  echo "Samples samengevoegd: $(ls *_segments.txt | wc -l)"
  echo "Regels in ALL_segments.txt: $(wc -l < ALL_segments.txt)"
  ;;

# STAP 6: genome-wide plots per sample in eigen map
plots)
  cd $WORK
  RF=$(mktemp --suffix=.R); trap 'rm -f $RF' EXIT
  cat > $RF << 'PLOTR'
chrs <- c(paste0("chr", 1:19), "chrX", "chrY")
cols <- c("1"="red3", "2"="orange", "3"="grey40", "4"="skyblue3", "5"="blue", "6"="purple4")
files <- list.files("results", pattern = "\\.rds$", full.names = TRUE)
message("Aantal samples: ", length(files))
for (f in files) {
  nm  <- sub("\\.rds$", "", basename(f))
  dir <- file.path("results", "plots", nm)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  orig <- file.path("results", "plots", paste0(nm, ".pdf"))
  if (file.exists(orig)) file.rename(orig, file.path(dir, paste0(nm, "_original.pdf")))
  tryCatch({
    x  <- readRDS(f)
    d  <- as.data.frame(x$corrected); d$chr <- as.character(d$chr)
    sg <- as.data.frame(x$segs$segs); sg$chr <- as.character(sg$chr)
    len <- as.numeric(tapply(d$end, d$chr, max)[chrs])
    off <- c(0, cumsum(len))[seq_along(chrs)]; names(off) <- chrs
    pdf(file.path(dir, paste0(nm, "_genomewide.pdf")), width = 14, height = 5)
    plot(d$start + off[d$chr], d$copy, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.35),
         ylim = c(-2, 2), xaxt = "n", xlab = "", ylab = "log2 copy", main = nm)
    abline(v = c(off, sum(len)), col = "grey85")
    axis(1, at = off + len / 2, labels = sub("chr", "", chrs), tick = FALSE, cex.axis = 0.7)
    segments(sg$start + off[sg$chr], sg$median, sg$end + off[sg$chr], sg$median,
             col = cols[as.character(sg$state)], lwd = 3)
    abline(h = 0, lty = 2)
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE)
                           message("FOUT bij ", nm, ": ", conditionMessage(e)) })
}
message("Klaar")
PLOTR
  module purge
  module load R-bundle-Bioconductor/3.23-foss-2026.1-R-4.6.1
  export OPENBLAS_NUM_THREADS=1
  Rscript $RF
  ;;

*)
  echo "Onbekende stap: $STEP"; exit 1
  ;;
esac