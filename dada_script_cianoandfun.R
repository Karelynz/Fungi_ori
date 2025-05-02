library(dada2)
library(ggplot2)
library(phyloseq)
library(Biostrings)
library(httr)

# Definir URL del repositorio de GitHub
repo_url <- "https://raw.githubusercontent.com/Karelynz/Fungi_ori/005ce89ea42129ebd14113819702dc3a801a1c66/procesadas/"

# Definir nombres de archivos
files <- c("_ciano_R1.fastq.gz", "_ciano_R2.fastq.gz", "_hongos_R1.fastq.gz", "_hongos_R2.fastq.gz")

dir_path <- "./procesadas/"
dir.create(dir_path, showWarnings = FALSE)

# Descargar archivos
for (file in files) {
  download.file(paste0(repo_url, file), destfile = file.path(dir_path, file))
}

# Obtener archivos FASTQ
fnFs_ciano <- list.files(dir_path, pattern="_ciano_R1.fastq.gz", full.names=TRUE)
fnRs_ciano <- list.files(dir_path, pattern="_ciano_R2.fastq.gz", full.names=TRUE)
fnFs_hongos <- list.files(dir_path, pattern="_hongos_R1.fastq.gz", full.names=TRUE)
fnRs_hongos <- list.files(dir_path, pattern="_hongos_R2.fastq.gz", full.names=TRUE)

# Visualizar calidad de secuencias de 5 en 5
for (i in seq(1, length(fnFs_ciano), by=5)) {
  print(plotQualityProfile(fnFs_ciano[i:min(i+4, length(fnFs_ciano))]))
  print(plotQualityProfile(fnRs_ciano[i:min(i+4, length(fnRs_ciano))]))
}
for (i in seq(1, length(fnFs_hongos), by=5)) {
  print(plotQualityProfile(fnFs_hongos[i:min(i+4, length(fnFs_hongos))]))
  print(plotQualityProfile(fnRs_hongos[i:min(i+4, length(fnRs_hongos))]))
}

# Salida de archivos filtrados
filtFs_ciano <- file.path(dir_path, "filtered", basename(fnFs_ciano))
filtRs_ciano <- file.path(dir_path, "filtered", basename(fnRs_ciano))
filtFs_hongos <- file.path(dir_path, "filtered", basename(fnFs_hongos))
filtRs_hongos <- file.path(dir_path, "filtered", basename(fnRs_hongos))

# Filtrado y recorte
filterAndTrim(fnFs_ciano, filtFs_ciano, fnRs_ciano, filtRs_ciano, truncLen=c(240, 200), maxN=0, maxEE=c(2, 2), truncQ=2, rm.phix=TRUE, compress=TRUE)
filterAndTrim(fnFs_hongos, filtFs_hongos, fnRs_hongos, filtRs_hongos, truncLen=c(240, 200), maxN=0, maxEE=c(2, 2), truncQ=2, rm.phix=TRUE, compress=TRUE)

# Aprender errores
errF_ciano <- learnErrors(filtFs_ciano)
errR_ciano <- learnErrors(filtRs_ciano)
errF_hongos <- learnErrors(filtFs_hongos)
errR_hongos <- learnErrors(filtRs_hongos)

# Inferencia de secuencias
dadaFs_ciano <- dada(filtFs_ciano, err=errF_ciano)
dadaRs_ciano <- dada(filtRs_ciano, err=errR_ciano)
dadaFs_hongos <- dada(filtFs_hongos, err=errF_hongos)
dadaRs_hongos <- dada(filtRs_hongos, err=errR_hongos)

# Unir secuencias forward y reverse
mergers_ciano <- mergePairs(dadaFs_ciano, filtFs_ciano, dadaRs_ciano, filtRs_ciano)
mergers_hongos <- mergePairs(dadaFs_hongos, filtFs_hongos, dadaRs_hongos, filtRs_hongos)

# Construcción de la tabla de ASVs
seqtab_ciano <- makeSequenceTable(mergers_ciano)
seqtab_hongos <- makeSequenceTable(mergers_hongos)

# Eliminación de quimeras
seqtab_ciano_nochim <- removeBimeraDenovo(seqtab_ciano, method="consensus")
seqtab_hongos_nochim <- removeBimeraDenovo(seqtab_hongos, method="consensus")

# Cargar base de datos para la asignación taxonómica
# Ajustar la base de datos según la región del gen objetivo
cianoseq_db <- "path/to/cianoseq_database.fa.gz"
unite_db <- "path/to/sh_general_release_dynamic.fasta"

# Asignar taxonomía
taxonomy_ciano <- assignTaxonomy(seqtab_ciano_nochim, cianoseq_db)
taxonomy_hongos <- assignTaxonomy(seqtab_hongos_nochim, unite_db)

# Guardar resultados
write.csv(taxonomy_ciano, file="taxonomy_ciano.csv")
write.csv(taxonomy_hongos, file="taxonomy_hongos.csv")

# Imprimir resumen
print("Proceso finalizado. Archivos de taxonomía guardados.")


# Evaluación de la comunidad MOCK
mock_species <- c("Saccharomyces cerevisiae", "Cryptococcus neoformans")
mock_abundance <- c(12, 2)  # Valores esperados según el reporte

mock_tax_hongos <- taxonomy_hongos[rownames(taxonomy_hongos) %in% mock_species, ]
mock_observed_abundance <- colSums(seqtab_hongos_nochim[rownames(seqtab_hongos_nochim) %in% rownames(mock_tax_hongos), ])

mock_comparison <- data.frame(
  Species = mock_species,
  Expected_Abundance = mock_abundance,
  Observed_Abundance = mock_observed_abundance
)
write.csv(mock_comparison, file="mock_evaluation_hongos.csv")

# 🔹 FILTRADO POR ABUNDANCIA DESPUÉS DE LA EVALUACIÓN MOCK 🔹
# 1. Eliminar singletons y doubletons
seqtab_ciano_nochim <- seqtab_ciano_nochim[, colSums(seqtab_ciano_nochim) > 2]
seqtab_hongos_nochim <- seqtab_hongos_nochim[, colSums(seqtab_hongos_nochim) > 2]

# 2. Corte basado en Bokulich et al. 2012 (0.005% del promedio de secuencias por muestra)
mean_reads_ciano <- mean(rowSums(seqtab_ciano_nochim))
min_abundance_ciano <- round(0.00005 * mean_reads_ciano)

mean_reads_hongos <- mean(rowSums(seqtab_hongos_nochim))
min_abundance_hongos <- round(0.00005 * mean_reads_hongos)

seqtab_ciano_nochim <- seqtab_ciano_nochim[, colSums(seqtab_ciano_nochim) >= min_abundance_ciano]
seqtab_hongos_nochim <- seqtab_hongos_nochim[, colSums(seqtab_hongos_nochim) >= min_abundance_hongos]

# 3. Generar una tabla filtrada para visualización (>10% de abundancia total)
threshold_ciano <- 0.1 * sum(seqtab_ciano_nochim)
seqtab_ciano_top <- seqtab_ciano_nochim[, colSums(seqtab_ciano_nochim) > threshold_ciano]

threshold_hongos <- 0.1 * sum(seqtab_hongos_nochim)
seqtab_hongos_top <- seqtab_hongos_nochim[, colSums(seqtab_hongos_nochim) > threshold_hongos]

# Guardar tablas filtradas
write.csv(seqtab_ciano_nochim, file="seqtab_ciano_filtered.csv")
write.csv(seqtab_hongos_nochim, file="seqtab_hongos_filtered.csv")
write.csv(seqtab_ciano_top, file="seqtab_ciano_top.csv")
write.csv(seqtab_hongos_top, file="seqtab_hongos_top.csv")

# 🔹 CREACIÓN DEL OBJETO PHYLOSEQ 🔹
# Crear objetos phyloseq para análisis en R
ps_ciano <- phyloseq(
  otu_table(seqtab_ciano_nochim, taxa_are_rows=FALSE),
  tax_table(as.matrix(taxonomy_ciano))
)

ps_hongos <- phyloseq(
  otu_table(seqtab_hongos_nochim, taxa_are_rows=FALSE),
  tax_table(as.matrix(taxonomy_hongos))
)

# Guardar objetos phyloseq
saveRDS(ps_ciano, file="phyloseq_ciano.rds")
saveRDS(ps_hongos, file="phyloseq_hongos.rds")

# Imprimir resumen
print("Proceso finalizado. Evaluación MOCK realizada, filtrado aplicado y objetos phyloseq creados.")
