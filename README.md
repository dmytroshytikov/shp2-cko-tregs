# SHP-2 Treg RNA-seq analysis

Analysis code for Figure 3 of "SHP-2 deficiency enhances the immunosuppressive
function of CD4+ regulatory T cells."

The comparison includes three stimulated WT and three stimulated Ptpn11-cKO
Treg samples. Positive log2 fold change indicates higher expression in cKO. The
workflow uses processed differential-expression statistics and sample-level
expression values.

## Run

    Rscript analysis.R

To specify an output directory:

    Rscript analysis.R /path/to/results

Input paths are resolved relative to `analysis.R`. All required R packages must
be installed before running the script.

## File contents

- `analysis.R`: performs GO analysis and creates numerical source tables and
  reproduces Figure 3A–E plots.
- `data/processed_expression.tsv.gz`: gene identifiers, Entrez mappings,
  differential-expression statistics, and six sample-expression columns. This
  dataset was derived from `data/original-data/DiffExpr.annotated.xls` after
  manual removal of the rows that have `NA` values in the `FDR` column.
- `data/samples.tsv`: sample identifiers, display labels, and genotypes.
- `data/original-data/DiffExpr.annotated.xls`: gene identifiers, Entrez
  mapping, normalized differential expression statistics (the raw read counts
  divided by the DESeq2 sample-specific size factors), six sample-expression
  columns, GO terms. This dataset includes all genes without removing invalid
  gene entries. Rows with `NA` values in the `FDR` column were manually removed
  to create `processed_expression.tsv.gz`.
- `sessionInfo.txt`: R and package versions.
- `SHA256SUMS.txt`: file checksums.

## Analysis settings

- Strict DEGs (Differentially Expressed Genes): adjusted P < 0.05 and absolute
  log2 fold change > 1.
- PCA: log2(expression + 1), the 1,000 most variable genes, centering without
  feature scaling.
- Volcano plot: adjusted P < 0.05 and absolute log2 fold change > 0.5 for
  coloring; -log10(adjusted P) capped at 60. Selected labels are independent of
  this display cutoff.
- Heatmap: Ptpn11 and the 39 highest-ranked eligible non-Ptpn11 protein-coding
  DEGs, excluding missing symbols and predicted Gm symbols.  Row-standardized
  log2 expression, Euclidean distance, and Ward D2 clustering are used. Color
  values are capped at +/-2.5.
- Directional GO over-representation: upregulated and downregulated strict DEGs
  analyzed separately with clusterProfiler::enrichGO.  BP, CC, and MF are
  tested with gene-set sizes of 10–500, BH correction, and adjusted P < 0.05.
  Figure 3D displays selected significant terms.
- Selected-gene panel: log2 fold changes for the specified genes; colors
  indicate adjusted P < 0.05 without a fold-change threshold.

## Analysis summary

- Tested genes: 11,919.
- Strict DEGs: 106, including 66 higher and 40 lower in cKO.
- Unique mapped Entrez identifiers supplied as the GO universe: 11,670.
- Effective GO annotation backgrounds: BP, 10,784; CC, 10,906; MF, 10,658.

## Outputs

- `results/tables/`: DEG lists, directional GO results, ontology background
  sizes, PCA coordinates, heatmap values and ordering, and selected-gene effect
  sizes.
- `results/figures/`: vector PDF and PNG files for panels A–E and an assembled
  analysis figure. The figure is further modified in Adobe Illustrator to
  improve visual presentation.
- `results/sessionInfo.txt`: the environment used for the analysis run.

## Software

R 4.5.3; clusterProfiler 4.18.4; org.Mm.eg.db 3.22.0; GO.db 3.22.0;
AnnotationDbi 1.72.0; dplyr 1.2.1; ggplot2 4.0.3; ggrepel 0.9.8;
patchwork 1.3.2; scales 1.4.0; tidyr 1.3.2.

The script requires org.Mm.eg.db and GO.db version 3.22.0.
Full environment details are provided in `sessionInfo.txt`.
