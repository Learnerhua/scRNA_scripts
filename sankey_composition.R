#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Composition sankey / alluvial plot (group x cell type).
#
# Companion helper for scRNA_composition_sccoda.py (step 8).  Kept as a
# standalone R script so the alluvial look stays faithful to the reference
# implementation (ggalluvial) instead of being re-implemented in Python.
#
# Usage:
#   Rscript sankey_composition.R <in_csv> <out_png> <out_pdf> <width> <height> \
#           [<group_order>]
#
# <in_csv>  : CSV whose first column is the group and whose remaining columns
#             are cell types; values are proportions (each row sums to 1).
#             **The column order IS the cell-type order** -- write the CSV with
#             the columns already in the order you want (the caller puts them in
#             global-abundance order, with an optional designated type last).
# <group_order>: optional comma-separated group order in LOGICAL order (e.g.
#                "control,treatment").  The strata are drawn bottom-to-top, so
#                the order is reversed internally: the FIRST entry ends up on
#                TOP.  (This mirrors the reference, which calls the plotter
#                with rev(group_order).)
#
# Requires: ggplot2, dplyr, tidyr, ggalluvial, RColorBrewer, scales
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(ggalluvial); library(RColorBrewer); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("usage: sankey_composition.R <in_csv> <out_png> <out_pdf> <width> <height> [group_order]")
}
in_csv  <- args[1]
out_png <- args[2]
out_pdf <- args[3]
fig_w   <- as.numeric(args[4])
fig_h   <- as.numeric(args[5])
group_order <- if (length(args) >= 6 && nzchar(args[6])) {
  strsplit(args[6], ",", fixed = TRUE)[[1]]
} else {
  character(0)
}

prop_data <- read.csv(in_csv, check.names = FALSE, stringsAsFactors = FALSE)
group_col <- colnames(prop_data)[1]

# -- long format -------------------------------------------------------------
plot_data <- prop_data %>%
  pivot_longer(cols = -all_of(group_col), names_to = "CellType", values_to = "Freq")

# Cell-type order = the CSV column order (see the usage note above).
celltype_levels <- colnames(prop_data)[-1]
if (length(group_order) > 0) {
  # Reverse for display: the first logical entry (control) goes on TOP.
  group_order <- rev(group_order)
  plot_data <- plot_data %>%
    filter(.data[[group_col]] %in% group_order) %>%
    mutate("{group_col}" := factor(.data[[group_col]], levels = group_order))
}
plot_data$CellType <- factor(plot_data$CellType, levels = rev(celltype_levels))

# -- palette: brewer Paired, extended by colorRampPalette when needed --------
# NB: faithful to the reference implementation -- colours are keyed by the
# ALPHABETICAL cell-type order (so a given cell type keeps its colour no matter
# how the abundance ordering changes), while the legend/stratum order follows
# celltype_levels above.
color_max <- brewer.pal.info["Paired", "maxcolors"]
nCluster  <- ncol(prop_data) - 1
base_colors <- brewer.pal(min(nCluster, color_max), "Paired")
palette_colors <- if (nCluster > color_max) {
  colorRampPalette(base_colors)(nCluster)
} else {
  base_colors
}
names(palette_colors) <- sort(colnames(prop_data)[-1])

# -- sankey ------------------------------------------------------------------
p <- ggplot(plot_data,
            aes(x = .data[[group_col]], stratum = CellType, alluvium = CellType,
                y = Freq, fill = CellType, label = CellType)) +
  scale_fill_manual(values = palette_colors, breaks = celltype_levels, name = NULL) +
  geom_flow(stat = "alluvium", lode.guidance = "leftright",
            color = "white", linewidth = 1.5) +
  geom_stratum(color = "transparent") +
  ylab("Frequency (%)") +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 12, face = "bold"),
    axis.text.y = element_text(color = "black", size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks.x = element_line(linewidth = 1),
    axis.ticks.length = unit(-4, "pt"),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

ggsave(out_png, p, width = fig_w, height = fig_h, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = fig_w, height = fig_h)
cat("sankey written:", out_png, "\n")
