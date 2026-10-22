###############################################################################
# CONTROL DE CALIDAD ATAC-seq: FRiP y enriquecimiento en TSS (hg38)
# Una muestra paired-end por ejecución.
###############################################################################

library(Rsamtools)
library(GenomicRanges)
library(rtracklayer)
library(Rsubread)
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

# Archivos de una misma réplica: ajustar los nombres según los archivos reales.
sample_id <- "SRR29829579"
bam_file <- "data/ATAC/SRR29829579.bam"
peaks_path <- "data/ATAC/SRR29829579_Control_rep1.bed"

###############################################################################
# 1. FRiP: Fraction of Fragments in Peaks
###############################################################################
# Definición operativa: fracción de pares primarios, no duplicados, con MAPQ >=30,
# ambos extremos mapeados, que solapan al menos un pico de la misma muestra.
# Numerador y denominador se obtienen con featureCounts y los mismos filtros.
# El denominador se calcula asignando fragmentos al conjunto de cromosomas
# utilizados en el análisis (no con countBam()/2).

peaks <- rtracklayer::import(peaks_path, format = "narrowPeak")
peaks <- GenomeInfoDb::keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- GenomicRanges::reduce(peaks, ignore.strand = TRUE)

bam_header <- Rsamtools::scanBamHeader(bam_file)[[1]]
chr_lengths <- bam_header$targets
common_chr <- intersect(names(chr_lengths), as.character(GenomeInfoDb::seqlevels(peaks)))
if (!length(common_chr)) stop("No coinciden los cromosomas del BAM y los peaks")
peaks <- GenomeInfoDb::keepSeqlevels(peaks, common_chr, pruning.mode = "coarse")

peaks_saf <- data.frame(
  GeneID = paste0("peak_", seq_along(peaks)),
  Chr = as.character(seqnames(peaks)),
  Start = start(peaks), End = end(peaks), Strand = "."
)

# SAF de cromosomas completos: misma selección cromosómica para denominador.
genome_saf <- data.frame(
  GeneID = paste0("chr_", seq_along(common_chr)),
  Chr = common_chr,
  Start = 1L,
  End = as.integer(chr_lengths[common_chr]),
  Strand = "."
)

count_fragments <- function(annotation) {
  Rsubread::featureCounts(
    files = bam_file,
    annot.ext = annotation,
    isGTFAnnotationFile = FALSE,
    isPairedEnd = TRUE,
    countReadPairs = TRUE,
    requireBothEndsMapped = TRUE,
    countMultiMappingReads = FALSE,
    primaryOnly = TRUE,
    ignoreDup = TRUE,
    minMQS = 30,
    allowMultiOverlap = TRUE,
    nthreads = 4
  )
}

# Cada fragmento puede asignarse a varios peaks; para FRiP usamos la cantidad
# de fragmentos 'Assigned' de la estadística, NO sumamos la matriz de counts.
get_assigned <- function(fc) {
  stat <- fc$stat
  if (!"Assigned" %in% stat$Status) stop("No se encontró 'Assigned' en featureCounts$stat")
  as.numeric(stat[stat$Status == "Assigned", 2])
}

fc_peaks <- count_fragments(peaks_saf)
fc_genome <- count_fragments(genome_saf)
n_fragments_in_peaks <- get_assigned(fc_peaks)
n_fragments_usable <- get_assigned(fc_genome)
if (n_fragments_usable == 0) stop("No se contaron fragmentos utilizables; revisar BAM y filtros")
frip <- n_fragments_in_peaks / n_fragments_usable
if (frip > 1) warning("FRiP > 1: revisar asignación, flags y definición de fragmento")

frip_summary <- data.frame(
  sample = sample_id,
  fragments_in_peaks = n_fragments_in_peaks,
  usable_fragments = n_fragments_usable,
  FRiP = frip
)

print(frip_summary)
write.csv(frip_summary, file.path("results", paste0(sample_id, "_FRiP.csv")), row.names = FALSE)

###############################################################################
# 2. Enriquecimiento de señal ATAC-seq alrededor del TSS
###############################################################################
# IMPORTANTE: el TSS enrichment de ATAC es una métrica de inserciones Tn5.
# TSSEscore() debe comprobarse en la versión instalada de ATACseqQC: su
# estructura de retorno y su compatibilidad con los BAM filtrados varían.

txs <- GenomicFeatures::transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
tss_upstream <- 1000
tss_downstream <- 1000
tss_end_size <- 100
tss_window <- 100

# Alinear nombres de cromosomas entre la anotación y el BAM.
common_tss_chr <- intersect(GenomeInfoDb::seqlevels(txs), names(chr_lengths))
if (!length(common_tss_chr)) stop("No coinciden los cromosomas del BAM y TxDb hg38")
txs <- GenomeInfoDb::keepSeqlevels(txs, common_tss_chr, pruning.mode = "coarse")

tss_regions <- GenomicRanges::promoters(
  txs, upstream = tss_upstream, downstream = tss_downstream
)
tss_regions_reduced <- GenomicRanges::reduce(tss_regions)

tss_alignments <- ATACseqQC::readBamFile(
  bamFile = bam_file,
  tag = character(0),
  which = tss_regions_reduced,
  asMates = FALSE,
  bigFile = TRUE
)

tss_result <- ATACseqQC::TSSEscore(
  obj = tss_alignments,
  txs = txs,
  upstream = tss_upstream,
  downstream = tss_downstream,
  endSize = tss_end_size,
  width = tss_window,
  step = tss_window,
  bigFile = TRUE
)

# Inspeccionar la salida ANTES de asumir el formato del perfil.
print(names(tss_result))
str(tss_result, max.level = 1)

if (is.null(tss_result$TSSEscore)) {
  warning("No se encontró $TSSEscore. Revisar la estructura de tss_result")
} else {
  tss_score <- tss_result$TSSEscore
  print(tss_score)
  write.csv(data.frame(sample = sample_id, TSS_enrichment = tss_score),
            file.path("results", paste0(sample_id, "_TSS_score.csv")), row.names = FALSE)
}

# Graficar SOLO si $values representa un vector numérico unidimensional.
# Verificar en la documentación de la versión instalada que sus puntos
# corresponden a la ventana y al paso configurados antes de interpretar el eje X.
if (!is.null(tss_result$values) && is.numeric(tss_result$values) &&
    is.null(dim(tss_result$values)) && length(tss_result$values) > 1) {

  tss_profile <- as.numeric(tss_result$values)
  tss_positions <- seq(-tss_upstream, tss_downstream,
                       length.out = length(tss_profile))
  tss_plot_df <- data.frame(position = tss_positions,
                            enrichment = tss_profile)
  p_tss <- ggplot(tss_plot_df, aes(x = position, y = enrichment)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "ATAC-seq enrichment around TSS",
         subtitle = paste0(sample_id, " | TSS score = ",
                           if (exists("tss_score")) round(tss_score, 2) else "NA"),
         x = "Approximate position relative to TSS (bp)",
         y = "Normalized ATAC-seq signal") +
    theme_classic(base_size = 13)
  print(p_tss)
  ggsave(file.path("plots", paste0(sample_id, "_TSS_enrichment.pdf")),
         plot = p_tss, width = 7, height = 5)
} else {
  message("No se generó el gráfico: comprobar la estructura de tss_result$values")
}
