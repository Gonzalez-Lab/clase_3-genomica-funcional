#BiocManager::install("DiffBind")
library(DiffBind)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(ggplot2)
library(dplyr)

# Importar metadata
samples_ChIP <- read.csv2("data/sample_sheet_ChIP.csv") 
# samples ATAC no tiene columna bamControl
# samples ChIP tiene bamControl (Input)

# Crear objeto DiffBind
db_ChIP <- dba(sampleSheet = samples_ChIP)

# Cuantificar reads/fragments en regiones comunes (paso LENTO)
db_ChIP <- dba.count(
  db_ChIP,
  minOverlap = 2,  # picos presentes en al menos 2 de las 6 muestras
  summits = 250    # ventanas de ~501 pb centradas en las cumbres de los picos
)

# Guardar objeto cuantificado para no repetir el conteo
saveRDS(db_ChIP, "results/db_ChIP_counted.rds")

# Definir contraste
db_ChIP <- dba.contrast(
  db_ChIP,
  design = ~ Condition,
  contrast = c("Condition", "SAHA", "Control")
)

# Análisis diferencial
db_ChIP <- dba.analyze(db_ChIP, method = DBA_DESEQ2)

# Guardar objeto analizado para reutilizarlo en clase
saveRDS(db_ChIP, "results/db_ChIP.rds")
# db_ChIP <- readRDS("results/db_ChIP.rds")

db_ChIP <- readRDS("results/db_ChIP.rds")
#db_ChIP <- readRDS("results/db_ChIP.rds")
#------------------------------

# Heatmap de correlaciones entre las muestras
plot(db_ChIP)

# PCA plot
dba.plotPCA(db_ChIP, attributes = DBA_CONDITION, label = DBA_ID)

# Volcano plot
dba.plotVolcano(db_ChIP)

# Extraer regiones diferenciales (FDR < 0.05)
differential_peaks_ChIP <- dba.report( db_ChIP,
                                       method = DBA_DESEQ2,
                                       th = 0.05 )

saveRDS(differential_peaks_ChIP, "results/differential_peaks_ChIP.rds")

# -----------------------------------------------------------------------------
# Anotar peaks con ChIPseeker
# -----------------------------------------------------------------------------
# Definición de región promotora alrededor del TSS
promoter_upstream   <- -3000
promoter_downstream <-  3000

diff_peak_annot_ChIP <- annotatePeak(
  peak = differential_peaks_ChIP,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  tssRegion = c(promoter_upstream, promoter_downstream),
  annoDb = "org.Hs.eg.db"
)

diff_peak_annot_ChIP

# Convertir resultado a data frame: conservar todos los picos con FDR < 0.05
diff_peak_annot_df_ChIP <- as.data.frame(diff_peak_annot_ChIP)
head(diff_peak_annot_df_ChIP)

# Verificar nombres de columnas disponibles
colnames(diff_peak_annot_df_ChIP)
stopifnot(all(c("FDR", "Fold", "annotation", "distanceToTSS") %in%
                colnames(diff_peak_annot_df_ChIP)))

# Guardar tabla completa de picos diferenciales anotados
write.csv(
  diff_peak_annot_df_ChIP,
  "results/differential_peaks_all_annotated_ChIP.csv",
  row.names = FALSE
)

# Seleccionar picos con FDR < 0.05 y |Fold| >= 1
# En DiffBind, Fold se reporta en escala log2; comprobar la dirección del contraste.
diff_filtered_ChIP <- diff_peak_annot_df_ChIP |>
  dplyr::filter(!is.na(FDR), !is.na(Fold), FDR < 0.05, abs(Fold) >= 1)

# Picos con mayor accesibilidad en SAHA
UP_ChIP <- diff_filtered_ChIP |> dplyr::filter(Fold >= 1)

# Picos con menor accesibilidad en SAHA
DOWN_ChIP <- diff_filtered_ChIP |> dplyr::filter(Fold <= -1)

cat("Picos diferenciales (FDR < 0.05):", nrow(diff_peak_annot_df_ChIP), "\n")
cat("Picos filtrados (FDR < 0.05 y |log2FC| >= 1):", nrow(diff_filtered_ChIP), "\n")
cat("UP:", nrow(UP_ChIP), "| DOWN:", nrow(DOWN_ChIP), "\n")

# Simplificar categorías de anotación para los picos filtrados
diff_filtered_ChIP$annotation_simple <- dplyr::case_when(
  grepl("^Promoter", diff_filtered_ChIP$annotation) ~ "Promoter",
  grepl("^5' UTR", diff_filtered_ChIP$annotation) ~ "5' UTR",
  grepl("^3' UTR", diff_filtered_ChIP$annotation) ~ "3' UTR",
  grepl("^Exon", diff_filtered_ChIP$annotation) ~ "Exon",
  grepl("^Intron", diff_filtered_ChIP$annotation) ~ "Intron",
  grepl("^Downstream", diff_filtered_ChIP$annotation) ~ "Downstream",
  grepl("Distal Intergenic", diff_filtered_ChIP$annotation) ~ "Distal intergenic",
  TRUE ~ "Other"
)

# Guardar tabla filtrada y anotada
write.csv(
  diff_filtered_ChIP,
  "results/differential peaks annotated _ChIP.csv",
  row.names = FALSE
)

# Resumen de anotación genómica de los picos filtrados
diff_annot_summary_ChIP <- diff_filtered_ChIP |>
  dplyr::count(annotation_simple, name = "n_peaks") |>
  dplyr::mutate(percent = n_peaks / sum(n_peaks) * 100) |>
  dplyr::arrange(dplyr::desc(n_peaks))

diff_annot_summary_ChIP

# Guardar summary
write.csv(
  diff_annot_summary_ChIP,
  "results/differential peak annotation summary _ChIP.csv",
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Gráfico de anotación genómica
# -----------------------------------------------------------------------------
p_annotation <- ggplot(
  diff_annot_summary_ChIP,
  aes(x = reorder(annotation_simple, percent), y = percent)
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Genomic annotation of differential peaks",
    subtitle = paste0(
      format(nrow(diff_filtered_ChIP), big.mark = ","),
      " differential peaks (FDR < 0.05, |log2FC| >= 1)"
    ),
    x = NULL,
    y = "Peaks (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

p_annotation

# Guardar plot
ggsave(
  "plots/differential peaks annotation _ChIP.pdf",
  plot = p_annotation,
  width = 7,
  height = 5
)

# -----------------------------------------------------------------------------
# Gráfico de distancia al TSS más cercano
# -----------------------------------------------------------------------------
distance_tss_plot_limit <- 100000 # 100 kb

p_distance_tss <- ggplot(
  diff_filtered_ChIP |>
    dplyr::filter(abs(distanceToTSS) <= distance_tss_plot_limit),
  aes(x = distanceToTSS)
) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Distance of differential peaks to nearest TSS",
    subtitle = "FDR < 0.05, |log2FC| >= 1",
    x = "Distance to nearest TSS (bp)",
    y = "Number of peaks"
  ) +
  theme_classic(base_size = 13)

p_distance_tss

# Guardar plot
ggsave(
  filename = "plots/diff peaks distance to TSS _ChIP.pdf",
  plot = p_distance_tss,
  width = 7,
  height = 5
)

###############################################################################
# Seleccionar picos diferenciales en promotores
###############################################################################

diff_promoter_peaks_ChIP <- diff_filtered_ChIP |>
  dplyr::filter(annotation_simple == "Promoter")

# Conservar información estadística
diff_promoter_peaks_df_ChIP <- diff_promoter_peaks_ChIP |>
  dplyr::select(
    seqnames,
    start,
    end,
    annotation,
    distanceToTSS,
    ENSEMBL,
    SYMBOL,
    GENENAME,
    Fold,
    FDR
  ) |>
  dplyr::filter(
    !is.na(ENSEMBL),
    !is.na(SYMBOL),
    SYMBOL != ""
  ) |>
  dplyr::arrange(Fold)

head(diff_promoter_peaks_df_ChIP)

# Picos en promotores con mayor accesibilidad en SAHA
promoter_up_ChIP <- diff_promoter_peaks_df_ChIP |>
  dplyr::filter(Fold >= 1)

# Picos en promotores con menor accesibilidad en SAHA
promoter_down_ChIP <- diff_promoter_peaks_df_ChIP |>
  dplyr::filter(Fold <= -1)

# Número de picos y genes únicos
cat("Picos en promotores UP:", nrow(promoter_up_ChIP), "\n")
cat("Picos en promotores DOWN:", nrow(promoter_down_ChIP), "\n")
cat("Genes únicos UP:", dplyr::n_distinct(promoter_up_ChIP$ENSEMBL), "\n")
cat("Genes únicos DOWN:", dplyr::n_distinct(promoter_down_ChIP$ENSEMBL), "\n")

# Guardar asociaciones peak-gen
write.csv(
  diff_promoter_peaks_df_ChIP,
  "results/differential promoter peaks _ChIP.csv",
  row.names = FALSE
)

# Guardar promotores UP y DOWN por separado
write.csv(promoter_up_ChIP, "results/promoter UP _ChIP.csv", row.names = FALSE)
write.csv(promoter_down_ChIP, "results/promoter DOWN _ChIP.csv", row.names = FALSE)
