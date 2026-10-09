# Readmapping en BAM-verwerking
# Werkmap instellen
setwd("C:/Users/user1/Desktop/Data_RA_raw/Data_RA_raw/")
getwd()
# Packages laden voor readmapping en BAM-verwerking
library(Rsubread)
library(Rsamtools)

# Softwareversies controleren
R.version.string
packageVersion("Rsubread")
packageVersion("Rsamtools")

# Referentiegenoom GRCh38 indexeren
# De index wordt gebruikt om RNA-seq-reads te aligneren.
buildindex(
  basename = "RA_project",
  reference = "GCF_000001405.26_GRCh38_genomic.fna",
  memory = 4000,
  indexSplit = TRUE)

# Paired-end FASTQ-reads mappen tegen GRCh38
# Vier RA-samples en vier controles, telkens met twee FASTQ-bestanden.
# De alignments worden opgeslagen als BAM-bestanden.

# RA-samples
align.RA1 <- align(
  index = "RA_project",
  readfile1 = "SRR4785979_1_subset40k.fastq",
  readfile2 = "SRR4785979_2_subset40k.fastq",
  output_file = "RA1.BAM")

align.RA2 <- align(
  index = "RA_project",
  readfile1 = "SRR4785980_1_subset40k.fastq",
  readfile2 = "SRR4785980_2_subset40k.fastq",
  output_file = "RA2.BAM")

align.RA3 <- align(
  index = "RA_project",
  readfile1 = "SRR4785986_1_subset40k.fastq",
  readfile2 = "SRR4785986_2_subset40k.fastq",
  output_file = "RA3.BAM")

align.RA4 <- align(
  index = "RA_project",
  readfile1 = "SRR4785988_1_subset40k.fastq",
  readfile2 = "SRR4785988_2_subset40k.fastq",
  output_file = "RA4.BAM")

# Controle-samples
align.ctrl1 <- align(
  index = "RA_project",
  readfile1 = "SRR4785819_1_subset40k.fastq",
  readfile2 = "SRR4785819_2_subset40k.fastq",
  output_file = "ctrl1.BAM")

align.ctrl2 <- align(
  index = "RA_project",
  readfile1 = "SRR4785820_1_subset40k.fastq",
  readfile2 = "SRR4785820_2_subset40k.fastq",
  output_file = "ctrl2.BAM")

align.ctrl3 <- align(
  index = "RA_project",
  readfile1 = "SRR4785828_1_subset40k.fastq",
  readfile2 = "SRR4785828_2_subset40k.fastq",
  output_file = "ctrl3.BAM")

align.ctrl4 <- align(
  index = "RA_project",
  readfile1 = "SRR4785831_1_subset40k.fastq",
  readfile2 = "SRR4785831_2_subset40k.fastq",
  output_file = "ctrl4.BAM")

# BAM-bestanden sorteren op genomische positie
# Gesorteerde BAM-bestanden zijn nodig voor efficiënte verwerking.
samples <- c(
  "RA1", "RA2", "RA3", "RA4",
  "ctrl1", "ctrl2", "ctrl3", "ctrl4")

lapply(samples, function(s) {
  sortBam(
    file = paste0(s, ".BAM"),
    destination = paste0(s, ".sorted")  )})

# BAM-bestanden indexeren
lapply(samples, function(s) {
  indexBam(file = paste0(s, ".sorted.bam"))})

# Count matrix maken met featureCounts

#Package laden voor het tellen van RNA-seq-reads
library(Rsubread)

# BAM-bestanden van de vier RA-samples en vier controles selecteren
allsamples <- c(
  "RA1.BAM", "RA2.BAM", "RA3.BAM", "RA4.BAM",
  "ctrl1.BAM", "ctrl2.BAM", "ctrl3.BAM", "ctrl4.BAM")

# Reads per gen tellen met featureCounts
# genomic.gtf bevat de genannotaties van NCBI RefSeq (GRCh38.p14).
# isPairedEnd = TRUE geeft aan dat paired-end sequencing is gebruikt.
# useMetaFeatures = TRUE groepeert exons op basis van gene_id.
count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE)

# Countgegevens uit het featureCounts-resultaat halen
str(count_matrix)
counts <- count_matrix$counts
head(counts)

# Kolomnamen aanpassen aan de bijbehorende samples
colnames(counts) <- c(
  "RA1", "RA2", "RA3", "RA4",
  "ctrl1", "ctrl2", "ctrl3", "ctrl4")

# Count matrix opslaan als CSV-bestand
# Rijen = genen; kolommen = samples; waarden = read counts.
write.csv(counts, "RA_countmatrix.csv")

# Differentiële genexpressieanalyse
# Packages laden
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)

# Count matrix inlezen
counts <- read.table(
  "count_matrix_RA.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE)

counts <- as.matrix(counts)

# Controleer de afmetingen, sample-ID's en gen-ID's.
dim(counts)
colnames(counts)
head(rownames(counts), 10)

# Samples indelen in RA en controle
# De SRA-run-ID's worden gekoppeld aan hun onderzoeksgroep.
groepen <- c(
  "SRR4785819" = "control",
  "SRR4785820" = "control",
  "SRR4785828" = "control",
  "SRR4785831" = "control",
  "SRR4785979" = "Reuma",
  "SRR4785980" = "Reuma",
  "SRR4785986" = "Reuma",
  "SRR4785988" = "Reuma")

# Controleer of alle sample-ID's correct overeenkomen.
stopifnot(setequal(colnames(counts), names(groepen)))

# DESeq2 vereist niet-negatieve, gehele read counts.
stopifnot(
  is.numeric(counts),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == round(counts)))

# Maak de sample-informatie voor DESeq2.
# 'control' wordt ingesteld als referentiegroep.
treatment_table <- data.frame(
  treatment = factor(
    groepen[colnames(counts)],
    levels = c("control", "Reuma") ),
  row.names = colnames(counts))
treatment_table
# DESeq2-analyse uitvoeren
# Het design vergelijkt genexpressie tussen RA en controles.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = treatment_table,
  design = ~ treatment)

# Normalisatie, dispersieschatting en statistische toetsing.
dds <- DESeq(dds)
# Positieve log2FoldChange betekent hogere expressie bij RA.
resultaten <- results(
  dds,
  contrast = c("treatment", "Reuma", "control")
)
resultsNames(dds)
# Significante genen identificeren
# padj < 0.05: Benjamini-Hochberg-gecorrigeerde p-waarde.
# |log2FoldChange| > 1: minstens een tweevoudig expressieverschil.

cat(
  "Totaal significant:",
  sum(resultaten$padj < 0.05, na.rm = TRUE),
  "\n")

cat(
  "Hoger bij RA:",
  sum(
    resultaten$padj < 0.05 &
      resultaten$log2FoldChange > 1,
    na.rm = TRUE  ),  "\n")

cat(
  "Lager bij RA:",
  sum(
    resultaten$padj < 0.05 &
      resultaten$log2FoldChange < -1,
    na.rm = TRUE  ),  "\n")

# Geselecteerde genen bekijken
# Controleer de expressieverschillen van biologisch
# relevante genen uit de resultaten.
genen <- c(
  "BCL2A1", "ADAMDEC1", "IGHV3-53",
  "IGHV1-69", "IGHG4", "IGHV4-31", "IGKV2-28")

res_df <- as.data.frame(resultaten)

res_df[
  intersect(genen, rownames(res_df)),
  c("log2FoldChange", "padj"),
  drop = FALSE]


# Resultaten rangschikken
# Bekijk genen met de hoogste en laagste log2FoldChange
# en genen met de kleinste aangepaste p-waarde.
hoogste_fold_change <- resultaten[
  order(resultaten$log2FoldChange, decreasing = TRUE),]

laagste_fold_change <- resultaten[
  order(resultaten$log2FoldChange),]

laagste_p_waarde <- resultaten[
  order(resultaten$padj),]

# Volcano plot maken
# De x-as toont de log2FoldChange en de y-as -log10(padj).
# De grenzen zijn padj < 0.05 en |log2FoldChange| > 1.
volcano <- EnhancedVolcano(
  resultaten,
  lab = rownames(resultaten),
  selectLab = c(
    "BCL2A1",
    "ADAMDEC1",
    "IGHV3-53",
    "IGHV1-69"  ),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1.2,
  labSize = 3.5,
  drawConnectors = TRUE,
  widthConnectors = 0.4,
  max.overlaps = Inf,
  legendPosition = "none",
  title = "Volcano plot",
  subtitle = "RA versus controle")

# Bekijk de volcano plot.
print(volcano)

# Sla de figuur op als PNG voor de README.
ggsave(
  filename = "VolcanoplotRA.png",
  plot = volcano,
  width = 12,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white")


#DESeq2-resultaten opslaan
# RDS-bestanden bewaren de R-objecten voor latere analyses.
saveRDS(dds, "dds_txt.rds")
saveRDS(resultaten, "resultaten_txt.rds")

# Controleer of beide bestanden zijn opgeslagen.
file.exists("dds_txt.rds")
file.exists("resultaten_txt.rds")

# Gene Ontology (GO)-verrijkingsanalyse

# Packages laden
library(goseq)
library(ggplot2)

# Significante genen selecteren
stopifnot(setequal(rownames(counts), rownames(resultaten)))

geteste_genen <- rownames(resultaten)[!is.na(resultaten$padj)]

# 1 = significant (padj < 0.05), 0 = niet significant.
gene_vector <- as.integer(
  resultaten[geteste_genen, "padj"] < 0.05
)
names(gene_vector) <- geteste_genen

table(gene_vector)


# Correctie voor selectiebias
countbias <- rowSums(counts)[geteste_genen]
stopifnot(!anyNA(countbias))

pwf <- nullp(
  gene_vector,
  bias.data = unname(countbias))

# GO-verrijkingsanalyse uitvoeren
# goseq onderzoekt welke GO-termen oververtegenwoordigd
# zijn onder de significant differentieel geëxpresseerde genen.
# De genen worden herkend aan hun gensymbolen.
GO_results <- goseq(
  pwf,
  genome = "hg19",
  id = "geneSymbol")

# Corrigeer de p-waarden voor meervoudig toetsen (BH).
GO_results$padj <- p.adjust(
  GO_results$over_represented_pvalue,
  method = "BH")

# Selecteer GO-termen met padj < 0.05.
sig_GO <- GO_results[
  !is.na(GO_results$padj) &
    GO_results$padj < 0.05,
]

sig_GO <- sig_GO[order(sig_GO$padj), ]

# Controleer het aantal significante GO-termen.
cat("Aantal significante GO-termen:", nrow(sig_GO), "\n")
head(sig_GO, 10)

# GO-resultaten opslaan
write.csv(
  GO_results,
  "GO_results_RA_all.csv",
  row.names = FALSE)

write.csv(
  sig_GO,
  "GO_results_RA_significant.csv",
  row.names = FALSE)

# Dotplot van de 10 significantste GO-termen
if (nrow(sig_GO) > 0) {
  
  # Selecteer de 10 GO-termen met de laagste padj.
  topGO <- head(sig_GO, 10)
  
  # Percentage DE-genen binnen iedere GO-term.
  topGO$hitsPerc <- 100 *
    topGO$numDEInCat / topGO$numInCat
  
  # Maak lange GO-termen beter leesbaar in de figuur.
  topGO$term_plot <- topGO$term
  
  # Verkort uitsluitend de langste term voor de figuur.
  topGO$term_plot[
    grepl(
      "^adaptive immune response based on",
      topGO$term_plot
    )
  ] <- "adaptive immune response (somatic recombination)"
  
  # Verdeel lange labels over meerdere regels.
  topGO$term_plot <- vapply(
    topGO$term_plot,
    function(x) {
      paste(strwrap(x, width = 32), collapse = "\n")
    },
    character(1) )
  
  # Maak de GO-dotplot.
  # Kleur = aangepaste p-waarde.
  # Puntgrootte = aantal DE-genen binnen de GO-term.
  go_plot <- ggplot(
    data = topGO,
    aes(
      x = hitsPerc,
      y = reorder(term_plot, hitsPerc),
      size = numDEInCat,
      colour = padj) ) +
    geom_point(alpha = 0.9) +
    scale_x_continuous(
      breaks = seq(0, 100, 20),
      limits = c(0, 100) ) +
    scale_colour_gradient(
      low = "darkblue",
      high = "dodgerblue",
      trans = "log10" ) +
    labs(
      title = "Top 10 GO terms (RA vs Control)",
      x = "Hits (%)",
      y = "GO term",
      size = "Number of DE genes",
      colour = "Adjusted p-value" ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 15   ),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 11),
      legend.position = "right",
      plot.margin = margin(15, 15, 15, 15) )
  
  # Figuur weergeven en opslaan.
  print(go_plot)
  
  ggsave(
    filename = "GO_plot_RA_definitief.png",
    plot = go_plot,
    width = 12,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"  )} else {message("Geen GO-termen met padj < 0.05 gevonden.")}

#  KEGG-pathwayvisualisatie met Pathview

# Packages laden
library(pathview)
library(org.Hs.eg.db)

# Gensymbolen omzetten naar Entrez-ID's
# Pathview gebruikt Entrez-ID's om genen te koppelen
# aan de juiste posities in de KEGG-pathway.
res_df <- as.data.frame(resultaten)

entrez_ids <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(res_df),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first")

# Log2-fold-changes koppelen aan Entrez-ID's
gene_fc <- res_df$log2FoldChange
names(gene_fc) <- entrez_ids

gene_fc <- gene_fc[
  !is.na(names(gene_fc)) &
    is.finite(gene_fc)]

# B-cell receptor signaling pathway visualiseren
# hsa04662 is de KEGG-code voor deze humane pathway.
# Positieve log2FC = hogere expressie bij RA (rood).
# Negatieve log2FC = lagere expressie bij RA (groen).
pathview(
  gene.data = gene_fc,
  pathway.id = "hsa04662",
  species = "hsa",
  gene.idtype = "ENTREZID",
  limit = list(gene = 5))