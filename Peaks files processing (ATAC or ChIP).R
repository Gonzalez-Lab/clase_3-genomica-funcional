library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(ggplot2)
library(dplyr)

#dir.create('results', showWarnings = FALSE)
#dir.create('plots', showWarnings = FALSE)

# Importar peaks replicas por condicion (galaxy saca narrowPeak)

peaks.1 <- import("data/SRR29829579_Control_rep1.bed", format = "narrowPeak")

peaks.2 <- import("data/SRR29829614_Control_rep2.bed", format = "narrowPeak")

peaks.3 <- import("data/SRR29829603_Control_rep3.bed", format = "narrowPeak")

# Definir replicas
peak_list <- list(peaks.1, peaks.2, peaks.3)

# Conservar cromosomas estandar
peak_list <- lapply(peak_list, function(x) {
  
  x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = "coarse" )
  reduce(x, ignore.strand = TRUE) #fusiona peaks solapados dentro de cada réplica
})

# Construir regiones candidatas, reduce une los picos solapados entre replicas
consensus_candidates <- reduce( do.call(c, peak_list),ignore.strand = TRUE )

# Contar replicas con evidencia de peak
support <- rowSums( sapply(peak_list, function(x) {
    countOverlaps(consensus_candidates, x, ignore.strand = TRUE) > 0})  )

# Criterio de reproducibilidad
min_replicates <- 2

# PICOS CONSENSO ENTRE TODAS LAS REPLICAS
consensus_peaks <- consensus_candidates[ support >= min_replicates ] 

length(consensus_peaks)

# Resumen descriptivo
peaks_summary <- data.frame(sample = 'ATAC-control', #poner el nombre de la condicion
                            n_peaks = length(consensus_peaks),
                            median_peak_width = median(width(consensus_peaks)),
                            mean_peak_width = mean(width(consensus_peaks)),
                            min_peak_width = min(width(consensus_peaks)),
                            max_peak_width = max(width(consensus_peaks))   )

print(peaks_summary)

write.csv(peaks_summary, 'results/ATAC control peaks summary.csv')

##############################################################################
# ANOTACION DE LOS PICOS
##############################################################################

# Definicion de region promotora alrededor del TSS
promoter_upstream   <- -3000
promoter_downstream <-  3000

# -----------------------------------------------------------------------------
# Anotar peaks con ChIPseeker
# -----------------------------------------------------------------------------
peak_annotation <- annotatePeak( peak = consensus_peaks,
                                 TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                 tssRegion = c(  promoter_upstream,
                                                 promoter_downstream ),
                                 annoDb = "org.Hs.eg.db" )

peak_annotation

# Convertir resultado a data frame
peak_annotation_df <- as.data.frame( peak_annotation )

head(peak_annotation_df)

# Simplificar categorias de anotacion
peak_annotation_df$annotation_simple <- dplyr::case_when(
  
  grepl("^Promoter", peak_annotation_df$annotation) ~ "Promoter",
  grepl("^5' UTR", peak_annotation_df$annotation) ~ "5' UTR",
  grepl("^3' UTR", peak_annotation_df$annotation) ~ "3' UTR",
  grepl("^Exon", peak_annotation_df$annotation) ~ "Exon",
  grepl("^Intron", peak_annotation_df$annotation) ~ "Intron",
  grepl("^Downstream", peak_annotation_df$annotation) ~ "Downstream",
  grepl("Distal Intergenic", peak_annotation_df$annotation) ~ "Distal intergenic",
  TRUE ~ "Other"
)

# Guardar tabla completa
write.csv( peak_annotation_df,"results/ATAC control consensus peaks annotated.csv")

annotation_summary <- peak_annotation_df |>
  dplyr::count( annotation_simple,  name = "n_peaks" ) |>
  dplyr::mutate( percent = n_peaks / sum(n_peaks) * 100 ) |>
  dplyr::arrange( dplyr::desc(n_peaks)  )

annotation_summary

#guardar summary
write.csv(annotation_summary, "results/ATAC_consensus_peak_annotation_summary.csv")

# -----------------------------------------------------------------------------
# Grafico de anotacion genomica
# -----------------------------------------------------------------------------

p_annotation <- ggplot(
  annotation_summary,
  aes( x = reorder(annotation_simple, percent), y = percent )
) + geom_col() + coord_flip() +
  labs(
    title = "Genomic annotation of reproducible ATAC-seq peaks",
    subtitle = paste0(
      format(length(consensus_peaks), big.mark = ",")," consensus peaks"),
    x = NULL,
    y = "Peaks (%)"
  ) +
  theme_classic(base_size = 13) +
  theme( plot.title = element_text(face = "bold") )

p_annotation

#guardar el plot
ggsave( "plots/ATAC control peaks annotation.pdf",
        plot = p_annotation,
        width = 7,  height = 5 )

# -----------------------------------------------------------------------------
# Grafico de distancia al TSS mas cercano
# -----------------------------------------------------------------------------

distance_tss_plot_limit <- 100000 # 100 kb

p_distance_tss <- ggplot( peak_annotation_df |>  dplyr::filter(
  abs(distanceToTSS) <= distance_tss_plot_limit  ),
  aes(x = distanceToTSS) ) +
  geom_histogram( bins = 100  ) +
  geom_vline(  xintercept = 0, linetype = "dashed" ) +
  labs(
    title = "Distance of reproducible ATAC-seq peaks to nearest TSS",
    x = "Distance to nearest TSS (bp)",
    y = "Number of peaks"
  ) +   theme_classic(base_size = 13)

p_distance_tss

# guardar plot
ggsave( filename = "plots/ATAC control distance to TSS.pdf",
        plot = p_distance_tss,
        width = 7, height = 5 )

#-------------------------------------------------------------------------------
# Seleccionar consensus peaks localizados en promotores
#------------------------------------------------------------------------------

promoter_peaks <- peak_annotation_df |> dplyr::filter(annotation_simple == "Promoter" )

# Construir tabla peak-gen
promoter_peak_df <- promoter_peaks |>
  dplyr::select(
    seqnames,
    start,
    end,
    annotation,
    distanceToTSS,
    ENSEMBL,
    SYMBOL,
    GENENAME ) |>  dplyr::filter(!is.na(ENSEMBL),
                                 !is.na(SYMBOL), 
                                 SYMBOL != "") |>
  #dplyr::distinct(geneId, SYMBOL, GENENAME) |>  
  dplyr::arrange(SYMBOL)

head(promoter_peak_df)

# Guardar asociaciones peak-gen
write.csv(promoter_peak_df,"results/ATAC control promoter_peak_gene_associations.csv")
