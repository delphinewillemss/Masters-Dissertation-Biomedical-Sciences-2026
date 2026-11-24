# CNV-achtige plots op basis van 500 kb bins: z = reads per bin / mediaan van het genoom
library(dplyr)
library(ggplot2)

base_dir   <- "/scratch/gent/vo/002/gvo00207/vsc52059"
input_file <- file.path(base_dir, "bin_counts", "bin_counts_all.tsv")
plot_dir   <- file.path(base_dir, "bin_counts", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# instellingen (nog af te stemmen met je mentor)
median_on_autosomes <- TRUE   # TRUE: mediaan over chr1-19; FALSE: over alle bins
loss_threshold <- 0.75        # z onder deze waarde = Loss
gain_threshold <- 1.25        # z boven deze waarde = Gain

# --- inlezen ---
d <- read.delim(input_file, stringsAsFactors = FALSE)
d$Start <- as.numeric(d$Start)
d$End   <- as.numeric(d$End)
d$Reads <- as.numeric(d$Reads) 
d$Chrom <- sub("^chr", "", as.character(d$Chrom))

chromosome_order <- c(as.character(1:19), "X", "Y")
chromosome_order <- chromosome_order[chromosome_order %in% unique(d$Chrom)]

d <- d[d$Chrom %in% chromosome_order, ]
d$Chrom <- factor(d$Chrom, levels = chromosome_order)

# --- mediaan en z per sample ---
d <- d %>%
  group_by(Sample) %>%
  mutate(
    Median = median(if (median_on_autosomes) Reads[Chrom %in% as.character(1:19)] else Reads),
    z = ifelse(Median > 0, Reads / Median, NA_real_)
  ) %>%
  ungroup()

# --- cumulatieve positie en chromosoomgrenzen ---
chr_info <- d %>%
  group_by(Chrom) %>%
  summarise(Length = as.numeric(max(End)), .groups = "drop") %>%
  arrange(Chrom) %>%
  mutate(
    Offset = cumsum(lag(Length, default = 0)),
    EndCum = Offset + Length,
    Midpoint = Offset + Length / 2
  )

d <- left_join(d, chr_info[, c("Chrom", "Offset")], by = "Chrom")
d$CumulativePosition <- (d$Start + d$End) / 2 + d$Offset

d$CNV <- cut(
  d$z,
  breaks = c(-Inf, loss_threshold, gain_threshold, Inf),
  labels = c("Loss", "Neutral", "Gain")
)

# --- tabel opslaan ---
write.table(
  d %>% select(Sample, Chrom, Start, End, Reads, Median, z, CNV),
  file.path(base_dir, "bin_counts", "bin_z_scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- plot per sample ---
for (s in unique(d$Sample)) {

  cat("Processing:", s, "\n")
  ds <- d[d$Sample == s & !is.na(d$z), ]

  p <- ggplot(ds, aes(CumulativePosition, z)) +
    geom_point(aes(color = CNV), size = 0.4, alpha = 0.6) +
    geom_hline(yintercept = 1, colour = "orange", linewidth = 0.4) +
    geom_vline(xintercept = chr_info$EndCum, colour = "grey80", linewidth = 0.3) +
    scale_color_manual(values = c(Loss = "red", Neutral = "black", Gain = "blue")) +
    scale_x_continuous(breaks = chr_info$Midpoint, labels = chr_info$Chrom) +
    coord_cartesian(ylim = c(0, 3)) +
    labs(
      title = paste("Read count ratio (500 kb bins) -", s),
      x = "Chromosome",
      y = "Reads / mediaan genoom"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust = 1)
    )

  ggsave(
    filename = file.path(plot_dir, paste0(s, "_bin_plot.png")),
    plot = p, width = 14, height = 5, dpi = 300
  )
}

cat("Done\n")