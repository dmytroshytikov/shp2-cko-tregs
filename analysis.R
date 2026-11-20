arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- arguments[startsWith(arguments, "--file=")]
if (length(script_argument) != 1) stop("Run this file with Rscript.")
root <- dirname(normalizePath(sub("^--file=", "", script_argument)))
output_arguments <- commandArgs(trailingOnly = TRUE)
if (length(output_arguments) > 1) stop("Usage: Rscript analysis.R [output_directory]")
result_dir <- if (length(output_arguments)) path.expand(output_arguments[1]) else file.path(root, "results")
output_dir <- file.path(result_dir, "figures")
tables_dir <- file.path(result_dir, "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("AnnotationDbi", "clusterProfiler", "org.Mm.eg.db", "GO.db",
                       "dplyr", "ggplot2", "ggrepel", "patchwork", "scales", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
for (package in c("org.Mm.eg.db", "GO.db")) {
  if (as.character(packageVersion(package)) != "3.22.0") {
    stop(package, " version 3.22.0 is required to preserve the manuscript's GO annotation.")
  }
}
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(tidyr)
})

set.seed(20260829)
data <- read.delim(gzfile(file.path(root, "data", "processed_expression.tsv.gz")),
                   check.names = FALSE, stringsAsFactors = FALSE)
samples <- read.delim(file.path(root, "data", "samples.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sample", "label", "genotype") %in% names(samples))) stop("Incomplete sample metadata.")
if (anyDuplicated(samples$sample) || anyDuplicated(samples$label)) stop("Sample identifiers must be unique.")
if (nrow(samples) != 6 || !identical(as.character(samples$genotype), c(rep("WT", 3), rep("cKO", 3)))) {
  stop("Expected three WT samples followed by three cKO samples.")
}
sample_columns <- samples$sample
sample_labels <- setNames(samples$label, samples$sample)
required_columns <- c("geneID", "ENSEMBL", "ENTREZID", "SYMBOL", "GENE_BIOTYPE",
                      "log2FoldChange", "pvalue", "padj", sample_columns)
if (!all(required_columns %in% names(data))) stop("Incomplete processed result table.")
if (anyDuplicated(data$geneID)) stop("Duplicate gene identifiers.")
if (any(!is.finite(data$pvalue) | !is.finite(data$padj) | !is.finite(data$log2FoldChange))) {
  stop("This frozen input should contain finite statistics for all tested genes.")
}
if (any(data$pvalue < 0 | data$pvalue > 1 | data$padj < 0 | data$padj > 1)) stop("Invalid P value.")
expression_values <- as.matrix(data[, sample_columns])
if (!is.numeric(expression_values) || any(!is.finite(expression_values) | expression_values < 0)) {
  stop("Invalid sample-expression values.")
}
data$ENTREZID <- as.character(data$ENTREZID)
data$ENTREZID[data$ENTREZID == ""] <- NA_character_

write_tsv <- function(table, filename) {
  write.table(table, file.path(tables_dir, filename), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "")
}
strict_hits <- data$padj < 0.05 & abs(data$log2FoldChange) > 1
strict_results <- data[strict_hits, ]
write_tsv(strict_results, "strict_DEGs.tsv")
universe_entrez <- unique(na.omit(data$ENTREZID))
gene_lists <- list(
  Up = unique(na.omit(data$ENTREZID[strict_hits & data$log2FoldChange > 0])),
  Down = unique(na.omit(data$ENTREZID[strict_hits & data$log2FoldChange < 0]))
)
run_go <- function(gene_ids, gene_set) {
  tables <- lapply(c("BP", "CC", "MF"), function(ontology) {
    enrichment <- suppressMessages(enrichGO(
      gene = gene_ids, universe = universe_entrez, OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID", ont = ontology, pAdjustMethod = "BH",
      pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10,
      maxGSSize = 500, readable = TRUE
    ))
    table <- as.data.frame(enrichment)
    if (!nrow(table)) return(NULL)
    ratio <- function(values) vapply(strsplit(as.character(values), "/", fixed = TRUE),
                                    function(parts) as.numeric(parts[1]) / as.numeric(parts[2]), numeric(1))
    table$RichFactor <- ratio(table$GeneRatio)
    table$FoldEnrichment <- table$RichFactor / ratio(table$BgRatio)
    table$GeneSet <- gene_set
    table$Ontology <- ontology
    table$AnnotationSource <- paste0("org.Mm.eg.db ", packageVersion("org.Mm.eg.db"))
    table
  })
  bind_rows(tables)
}
go_results <- bind_rows(Map(run_go, gene_lists, names(gene_lists))) %>%
  arrange(GeneSet, Ontology, p.adjust, pvalue)
if (!nrow(go_results)) stop("No GO results were returned.")
go <- go_results %>% filter(!is.na(p.adjust), p.adjust < 0.05)
write_tsv(go_results, "GO_all_results.tsv")
write_tsv(go, "GO_significant.tsv")
write_tsv(go_results %>%
  transmute(Ontology,
            input_mapped_universe = length(universe_entrez),
            annotated_background = as.integer(sub(".*/", "", BgRatio))) %>%
  distinct(), "GO_background_sizes.tsv")
write_tsv(data.frame(
  category = c("Tested genes", "Mapped universe", "Strict DEGs", "Up", "Down"),
  count = c(nrow(data), length(universe_entrez), sum(strict_hits),
            sum(strict_hits & data$log2FoldChange > 0),
            sum(strict_hits & data$log2FoldChange < 0))
), "analysis_counts.tsv")

save_plot <- function(plot, stem, width_mm, height_mm) {
  pdf_device <- function(filename, width, height, ...) {
    grDevices::pdf(filename, width = width, height = height,
                   family = "Helvetica", useDingbats = FALSE, ...)
  }
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
         width = width_mm, height = height_mm, units = "mm",
         device = pdf_device, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 300, bg = "white")
}

wt_color <- "#3B6FB6"
ko_color <- "#D95F4B"
up_color <- "#C85A37"
down_color <- "#3973AC"
ns_color <- "#B8B8B8"
heat_low <- "#3B6FB6"
heat_high <- "#D77A4A"
volcano_up_color <- "#D98D73"
volcano_down_color <- "#7E9EC4"
volcano_ns_color <- "#D7D7D7"

base_theme <- theme_classic(base_size = 7.6, base_family = "sans") +
  theme(
    plot.title = element_text(size = 8.7, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 6.4, color = "#4D4D4D", hjust = 0.5),
    axis.title = element_text(size = 7.4),
    axis.text = element_text(size = 6.6, color = "#222222"),
    legend.title = element_text(size = 6.3),
    legend.text = element_text(size = 6.2),
    legend.key.height = unit(2.8, "mm"),
    legend.key.width = unit(4.0, "mm"),
    strip.background = element_rect(fill = "#F3F3F3", color = "#C9C9C9", linewidth = 0.3),
    strip.text = element_text(size = 6.3, face = "bold"),
    plot.margin = margin(4, 5, 4, 6)
  )

expression_matrix <- as.matrix(data[, sample_columns])
storage.mode(expression_matrix) <- "numeric"
gene_labels <- ifelse(is.na(data$SYMBOL) | data$SYMBOL == "", data$ENSEMBL, data$SYMBOL)
rownames(expression_matrix) <- make.unique(gene_labels)
log_expression <- log2(expression_matrix + 1)

metadata <- as_tibble(samples) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "cKO")))

strict_de <- !is.na(data$padj) & data$padj < 0.05 & abs(data$log2FoldChange) > 1
display_de <- !is.na(data$padj) & data$padj < 0.05 & abs(data$log2FoldChange) > 0.5
data <- data %>%
  mutate(
    volcano_class = case_when(
      display_de & log2FoldChange < 0 ~ "Lower in cKO",
      display_de & log2FoldChange > 0 ~ "Higher in cKO",
      TRUE ~ "Below cutoff"
    ),
    volcano_class = factor(
      volcano_class,
      levels = c("Lower in cKO", "Below cutoff", "Higher in cKO")
    ),
    plot_padj = pmax(ifelse(is.na(padj), 1, padj), 1e-300),
    neg_log10_padj = pmin(-log10(plot_padj), 60)
  )

gene_variance <- apply(log_expression, 1, var, na.rm = TRUE)
top_variable <- order(gene_variance, decreasing = TRUE)[seq_len(min(1000, nrow(log_expression)))]
pca <- prcomp(t(log_expression[top_variable, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_variance <- 100 * summary(pca)$importance[2, 1:2]
pca_data <- as.data.frame(pca$x[, 1:2]) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(metadata, by = "sample") %>%
  mutate(
    label_x = PC1 + case_when(
      label == "WT2" ~ 0.8,
      label == "cKO2" ~ -0.5,
      TRUE ~ -0.9
    ),
    label_y = PC2 + case_when(
      label %in% c("WT1", "WT3", "cKO1", "cKO3") ~ -1.7,
      TRUE ~ 1.7
    ),
    label_hjust = ifelse(label == "WT2", 0, 1)
  )

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = genotype)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "#E1E1E1") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#E1E1E1") +
  geom_point(size = 2.7) +
  geom_text(
    aes(x = label_x, y = label_y, label = label, hjust = label_hjust),
    size = 2.15,
    show.legend = FALSE,
    vjust = 0.5
  ) +
  scale_color_manual(values = c(WT = wt_color, cKO = ko_color)) +
  labs(
    title = "Principal component analysis",
    x = sprintf("PC1 (%.1f%%)", pca_variance[1]),
    y = sprintf("PC2 (%.1f%%)", pca_variance[2]),
    color = NULL
  ) +
  base_theme +
  theme(legend.position = "top")

volcano_labels <- c("Ptpn11", "Foxp3", "Nt5e", "Cxcr4", "Il13", "Klrd1")
label_data <- data %>%
  filter(SYMBOL %in% volcano_labels) %>%
  mutate(plot_label = ifelse(SYMBOL == "Nt5e", "Nt5e/CD73", SYMBOL))

p_volcano <- ggplot(data, aes(log2FoldChange, neg_log10_padj, color = volcano_class)) +
  geom_point(size = 0.64) +
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed",
    linewidth = 0.3,
    color = "#737373"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.3,
    color = "#737373"
  ) +
  geom_point(
    data = label_data,
    shape = 21,
    fill = "white",
    color = "#333333",
    size = 1.65,
    stroke = 0.35
  ) +
  geom_label_repel(
    data = label_data,
    aes(label = plot_label),
    color = "#222222",
    fill = "white",
    label.size = NA,
    label.padding = unit(0.07, "lines"),
    size = 2.0,
    seed = 20260829,
    min.segment.length = 0,
    segment.size = 0.22,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.18,
    force = 2.5,
    force_pull = 0.5,
    max.time = 5,
    max.iter = 100000,
    ylim = c(1.8, 60)
  ) +
  scale_color_manual(values = c(
    "Lower in cKO" = volcano_down_color,
    "Below cutoff" = volcano_ns_color,
    "Higher in cKO" = volcano_up_color
  )) +
  coord_cartesian(xlim = c(-7.7, 7.7), ylim = c(0, 62), clip = "off") +
  labs(
    title = "Differential expression",
    x = "log2 fold change (cKO / WT)",
    y = expression(-log[10](adjusted~P)),
    color = NULL
  ) +
  base_theme +
  theme(
    legend.position = "top",
    legend.spacing.x = unit(0.7, "mm")
  )

eligible_heatmap_genes <- data %>%
  filter(
    !is.na(padj), padj < 0.05, abs(log2FoldChange) > 1,
    GENE_BIOTYPE == "protein_coding",
    !is.na(SYMBOL), SYMBOL != "", !grepl("^Gm[0-9]+", SYMBOL)
  ) %>%
  arrange(padj)

heatmap_candidates <- bind_rows(
  eligible_heatmap_genes %>%
    filter(SYMBOL != "Ptpn11") %>%
    slice_head(n = 39),
  eligible_heatmap_genes %>%
    filter(SYMBOL == "Ptpn11") %>%
    slice_head(n = 1)
) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  arrange(padj)

if (!"Ptpn11" %in% heatmap_candidates$SYMBOL) {
  stop("Ptpn11 was not available among strict protein-coding DEGs.")
}

top_matrix <- log2(as.matrix(heatmap_candidates[, sample_columns]) + 1)
rownames(top_matrix) <- heatmap_candidates$SYMBOL
top_z <- t(scale(t(top_matrix)))
top_z[is.na(top_z)] <- 0
row_order <- hclust(dist(top_z), method = "ward.D2")$order
column_order <- hclust(dist(t(top_z)), method = "ward.D2")$order
top_z <- top_z[row_order, column_order, drop = FALSE]

top_long <- as.data.frame(top_z) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
  mutate(
    sample = factor(sample, levels = colnames(top_z), labels = sample_labels[colnames(top_z)]),
    gene = factor(gene, levels = rev(rownames(top_z))),
    is_ptpn11 = as.character(gene) == "Ptpn11",
    z = pmax(pmin(z, 2.5), -2.5)
  )

p_heatmap <- ggplot(top_long, aes(sample, gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.16) +
  geom_tile(
    data = top_long %>% filter(is_ptpn11),
    fill = NA,
    color = "#222222",
    linewidth = 0.5
  ) +
  scale_fill_gradient2(
    low = heat_low,
    mid = "white",
    high = heat_high,
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = squish
  ) +
  guides(fill = guide_colorsteps(
    direction = "horizontal",
    title.position = "top",
    barwidth = unit(25, "mm"),
    barheight = unit(2.0, "mm"),
    even.steps = TRUE,
    show.limits = TRUE
  )) +
  labs(
    title = "Top protein-coding DEGs",
    x = NULL,
    y = NULL,
    fill = "Row z-score"
  ) +
  base_theme +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 5.8),
    axis.text.y = element_text(size = 4.9),
    panel.border = element_rect(color = "#777777", fill = NA, linewidth = 0.3),
    plot.margin = margin(4, 4, 3, 6)
  )

ora_terms <- tribble(
  ~GeneSet, ~ID, ~short_label,
  "Up", "GO:0034103", "Tissue remodeling",
  "Up", "GO:0048514", "Blood vessel morphogenesis",
  "Up", "GO:0006898", "Receptor-mediated endocytosis",
  "Up", "GO:0007160", "Cell-matrix adhesion",
  "Up", "GO:0009986", "Cell surface",
  "Down", "GO:0001906", "Cell killing",
  "Down", "GO:0045954", "Positive regulation of NK cytotoxicity",
  "Down", "GO:0001819", "Positive regulation of cytokine production",
  "Down", "GO:0070661", "Leukocyte proliferation",
  "Down", "GO:0005125", "Cytokine activity"
)

ora_selected <- ora_terms %>%
  left_join(go, by = c("GeneSet", "ID")) %>%
  mutate(
    signed_significance = ifelse(GeneSet == "Up", 1, -1) * -log10(p.adjust),
    direction_label = ifelse(GeneSet == "Up", "Higher in cKO", "Lower in cKO"),
    short_label = factor(short_label, levels = rev(short_label))
  )

if (any(is.na(ora_selected$p.adjust))) {
  stop("At least one selected GO term was not found in the enrichment table.")
}

p_go <- ggplot(ora_selected, aes(signed_significance, short_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "#555555") +
  geom_segment(
    aes(x = 0, xend = signed_significance, yend = short_label, color = direction_label),
    linewidth = 0.65
  ) +
  geom_point(aes(color = direction_label), size = 2.7) +
  scale_color_manual(values = c("Higher in cKO" = up_color, "Lower in cKO" = down_color)) +
  scale_x_continuous(labels = function(x) abs(x), expand = expansion(mult = c(0.03, 0.08))) +
  labs(
    title = "Directional gene-program shifts",
    x = expression("Lower in cKO  " %<-% "  " * -log[10](FDR) * "  " %->% "  Higher in cKO"),
    y = NULL,
    color = NULL
  ) +
  base_theme +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text.y = element_text(size = 5.8),
    plot.margin = margin(4, 7, 3, 5)
  )

focused_gene_design <- tribble(
  ~category, ~symbol, ~display, ~gene_order,
  "Treg identity", "Ptpn11", "Ptpn11", 1,
  "Treg identity", "Foxp3", "Foxp3", 2,
  "Treg identity", "Ctla4", "Ctla4", 3,
  "CD73 / GITR", "Nt5e", "Nt5e (CD73)", 1,
  "CD73 / GITR", "Entpd1", "Entpd1 (CD39)", 2,
  "CD73 / GITR", "Tnfrsf18", "Tnfrsf18 (GITR)", 3,
  "Suppressive cytokines", "Il10", "Il10", 1,
  "Suppressive cytokines", "Tgfb1", "Tgfb1", 2,
  "Suppressive cytokines", "Tgfb3", "Tgfb3", 3,
  "Suppressive cytokines", "Ebi3", "Ebi3 (IL-35 subunit)", 4,
  "Tr1-associated", "Ifng", "Ifng", 1,
  "Tr1-associated", "Il4", "Il4", 2,
  "Tr1-associated", "Lag3", "Lag3", 3,
  "Tr1-associated", "Itga2", "Itga2 (CD49b)", 4,
  "Tr1-associated", "Maf", "Maf", 5,
  "Tr1-associated", "Eomes", "Eomes", 6,
  "Tr1-associated", "Prdm1", "Prdm1", 7,
  "IL-2 / STAT5", "Il2ra", "Il2ra", 1,
  "IL-2 / STAT5", "Il2rb", "Il2rb", 2,
  "IL-2 / STAT5", "Il2rg", "Il2rg", 3,
  "IL-2 / STAT5", "Jak3", "Jak3", 4,
  "IL-2 / STAT5", "Stat5a", "Stat5a", 5,
  "IL-2 / STAT5", "Stat5b", "Stat5b", 6,
  "TGF-beta / SMAD", "Tgfbr1", "Tgfbr1", 1,
  "TGF-beta / SMAD", "Tgfbr2", "Tgfbr2", 2,
  "TGF-beta / SMAD", "Smad2", "Smad2", 3,
  "TGF-beta / SMAD", "Smad3", "Smad3", 4,
  "TGF-beta / SMAD", "Smad7", "Smad7", 5
)

category_levels <- c(
  "TGF-beta / SMAD", "IL-2 / STAT5", "Tr1-associated",
  "Suppressive cytokines", "CD73 / GITR", "Treg identity"
)

focused_gene_results <- focused_gene_design %>%
  left_join(
    data %>% select(SYMBOL, log2FoldChange, pvalue, padj),
    by = c("symbol" = "SYMBOL")
  ) %>%
  mutate(
    category = factor(category, levels = category_levels),
    effect_class = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Higher in cKO (FDR < 0.05)",
      !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Lower in cKO (FDR < 0.05)",
      TRUE ~ "Not significant"
    ),
    effect_class = factor(
      effect_class,
      levels = c("Lower in cKO (FDR < 0.05)", "Not significant", "Higher in cKO (FDR < 0.05)")
    )
  ) %>%
  arrange(category, gene_order) %>%
  mutate(plot_order = row_number())

if (any(is.na(focused_gene_results$log2FoldChange))) {
  stop("At least one focused gene was not available in the processed differential-expression table.")
}

write.table(
  focused_gene_results %>%
    select(plot_order, category, symbol, display, log2FoldChange, pvalue, padj, effect_class),
  file.path(tables_dir, "panelE_selected_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

focused_group_bands <- focused_gene_results %>%
  group_by(category) %>%
  summarise(
    xmin = min(plot_order) - 0.45,
    xmax = max(plot_order) + 0.45,
    xmid = mean(plot_order),
    .groups = "drop"
  ) %>%
  mutate(
    category_label = recode(
      as.character(category),
      "TGF-beta / SMAD" = "TGF-beta signaling",
      "IL-2 / STAT5" = "IL-2 / STAT5 signaling",
      "Tr1-associated" = "Tr1-associated",
      "Suppressive cytokines" = "Suppressive cytokines",
      "CD73 / GITR" = "CD73 / GITR",
      "Treg identity" = "Treg identity"
    )
  )

p_focused <- ggplot(focused_gene_results, aes(plot_order, log2FoldChange)) +
  geom_rect(
    data = focused_group_bands,
    aes(xmin = xmin, xmax = xmax, ymin = 1.18, ymax = 1.47),
    inherit.aes = FALSE,
    fill = "#F3F3F3",
    color = "#C9C9C9",
    linewidth = 0.28
  ) +
  geom_text(
    data = focused_group_bands,
    aes(x = xmid, y = 1.325, label = category_label),
    inherit.aes = FALSE,
    size = 1.65,
    fontface = "bold",
    color = "#333333"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.38, color = "#555555") +
  geom_segment(
    aes(xend = plot_order, y = 0, yend = log2FoldChange, color = effect_class),
    linewidth = 0.55
  ) +
  geom_point(aes(color = effect_class), size = 1.9) +
  scale_color_manual(values = c(
    "Lower in cKO (FDR < 0.05)" = down_color,
    "Not significant" = ns_color,
    "Higher in cKO (FDR < 0.05)" = up_color
  )) +
  scale_x_continuous(
    breaks = focused_gene_results$plot_order,
    labels = focused_gene_results$display,
    expand = expansion(add = c(0.65, 0.65))
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, 0.5),
    limits = c(-1.12, 1.49),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Selected immunoregulatory genes",
    subtitle = "Effect sizes are shown for all tested genes; gray indicates FDR at or above 0.05",
    x = NULL,
    y = "log2 fold change (cKO / WT)",
    color = NULL
  ) +
  base_theme +
  theme(
    legend.position = "top",
    legend.spacing.x = unit(0.8, "mm"),
    axis.text.y = element_text(size = 5.8),
    axis.text.x = element_text(size = 5.0, angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_line(linewidth = 0.25, color = "#777777"),
    plot.margin = margin(4, 6, 4, 7)
  )

final_figure <- wrap_plots(
  A = p_pca,
  B = p_volcano,
  C = p_heatmap,
  D = p_go,
  E = p_focused,
  design = "AAACCC\nAAACCC\nBBBCCC\nBBBCCC\nDDEEEE\nDDEEEE",
  widths = c(1, 1, 1, 1, 1, 1),
  heights = c(1, 1, 1, 1, 1, 1)
) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 9.5, family = "sans"),
      plot.tag.position = c(0.005, 0.995)
    )
  )


write_tsv(pca_data, "panelA_PCA_coordinates.tsv")
write_tsv(data %>% select(geneID, SYMBOL, log2FoldChange, pvalue, padj, volcano_class,
                          plot_padj, neg_log10_padj), "panelB_volcano_values.tsv")
write_tsv(heatmap_candidates %>% select(all_of(required_columns)), "panelC_heatmap_genes.tsv")
write.table(top_z, file.path(tables_dir, "panelC_unclipped_row_zscores.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write_tsv(data.frame(sample = colnames(top_z), label = sample_labels[colnames(top_z)],
                     order = seq_len(ncol(top_z))), "panelC_sample_order.tsv")
write_tsv(ora_selected, "panelD_selected_GO.tsv")
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
save_plot(p_pca, "Figure3A_PCA", 75, 70)
save_plot(p_volcano, "Figure3B_volcano", 95, 80)
save_plot(p_heatmap, "Figure3C_heatmap", 85, 150)
save_plot(p_go, "Figure3D_directional_GO", 105, 100)
save_plot(p_focused, "Figure3E_selected_genes", 183, 88)
save_plot(final_figure, "Figure3_ABCDE_analysis_preview", 183, 247)
message("Completed Figure 3 downstream analysis. Results: ", normalizePath(result_dir))

