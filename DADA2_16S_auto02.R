### DADA2 PIPELINE 2.0 WITH MAXIMUM MULTITHREADING ON WINDOWS

# Instalar librerC-as si es necesario
defaults <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("dada2", quietly = TRUE))
    BiocManager::install("dada2", version = "3.20")
}
defaults()

# Cargar librerC-as
library(dada2)
library(stats)
library(BiocParallel)
library(R.utils)
library(parallel)  # Para detecciC3n de nC:cleos

# Detectar nC:mero de nC:cleos disponibles
nc <- detectCores(logical = FALSE)
message("Usando ", nc, " nC:cleos fC-sicos para el procesamiento paralelo.")

# Configurar el backend paralelo (en Windows, usar SnowParam)
bpp <- SnowParam(workers = nc, type = "SOCK")
register(bpp, default = TRUE)

# Directorio de trabajo y archivos
dir <- "C:/Users/VladM/Documents/DEunir/1Cuatrim/seqyOmnics/Actividades/Actividad 3/seq_Act03WD_al03"
setwd(dir)

fnFs <- sort(list.files(dir, pattern = "_1.fastq.gz$", full.names = TRUE))
print(fnFs)
fnRs <- sort(list.files(dir, pattern = "_2.fastq.gz$", full.names = TRUE))
print(fnRs)

# Extraer nombres de muestras
auto_names <- function(files, suffix) {
  sub(suffix, "", basename(files))
}
baseF <- auto_names(fnFs, "_1.fastq.gz$")
baseR <- auto_names(fnRs, "_2.fastq.gz$")
paired <- intersect(baseF, baseR)
sample.names <- sapply(paired, function(x) paste(unlist(strsplit(x, "_"))[1:2], collapse = "_"))
print(sample.names)

# InspecciC3n de calidad (opcional)
plotQualityProfile(fnFs[1:min(12, length(fnFs))])
plotQualityProfile(fnRs[1:min(12, length(fnRs))])

# Crear directorio filtered si no existe
filt_path <- file.path(dir, "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)

# Generar nombres de archivos filtrados
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtRs) <- sample.names

# Filtrado y trimming con multithread = nC:mero de nC:cleos
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(200, 180),   # antes 240/200
                     maxN = 0, maxEE = c(2,2), truncQ = 2,
                     rm.phix = TRUE, compress = TRUE,
                     multithread = nc)

# Aprender errores en paralelo
errF <- learnErrors(filtFs, multithread = nc )
errR <- learnErrors(filtRs, multithread = nc)

# Visualizar errores (opcional)
plotErrors(errF, nominalQ = TRUE)

# DereplicaciC3n
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Inferencia DADA (pooled y no-pooled) en paralelo
dadaFs <- dada(derepFs, err = errF, multithread = nc)   # idem
dadaRs <- dada(derepRs, err = errR, multithread = nc)

dadaPFs <- dada(derepFs, err = errF, multithread = nc, pool = TRUE)
dadaPRs <- dada(derepRs, err = errR, multithread = nc, pool = TRUE)

# Merge paired reads
mergers   <- mergePairs(dadaFs, derepFs, dadaRs, derepRs  )
mergersP  <- mergePairs(dadaPFs, derepFs, dadaPRs, derepRs)

# Tablas de secuencias
seqtab    <- makeSequenceTable(mergers)
seqtabP   <- makeSequenceTable(mergersP)

# EliminaciC3n de quimeras en paralelo
seqtab.nochim  <- removeBimeraDenovo(seqtab, method = "consensus",  verbose = TRUE)
seqtab.nochimP <- removeBimeraDenovo(seqtabP, method = "consensus",  verbose = TRUE)

# Seguimiento de lecturas a travC)s del pipeline
getN <- function(x) sum(getUniques(x))
track  <- cbind(out,
                sapply(dadaFs, getN),
                sapply(dadaRs, getN),
                sapply(mergers, getN),
                rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names

trackP <- cbind(out,
                sapply(dadaPFs, getN),
                sapply(dadaPRs, getN),
                sapply(mergersP, getN),
                rowSums(seqtab.nochimP))
colnames(trackP) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(trackP) <- sample.names

trackP

# Guardar resultados
saveRDS(track, file = file.path(dir, "track_reads.rds"))
saveRDS(trackP, file = file.path(dir, "track_reads_pooled.rds"))
save.image(file = file.path(dir, "dada2_full_workspace.RData"))

#------------------------------------------------------------------------------------------------------
###ASSIGN TAXONOMY### (corriC3 la de la carpeta databases de todos)

# Ruta a la base de datos (descC!rgala previamente desde dada2 website)
silva_ref <- "C:/Users/VladM/Documents/DEunir/1Cuatrim/seqyOmnics/Actividades/Actividad 3/seq_Act03WD_al03/Bases_Datos/silva_nr99_v138.1_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab.nochimP, silva_ref, multithread = nc)
taxa <- addSpecies(taxa, "C:/Users/VladM/Documents/DEunir/1Cuatrim/seqyOmnics/Actividades/Actividad 3/seq_Act03WD_al03/Bases_Datos/silva_species_assignment_v138.1.fa.gz") # Opcional



# Guardar la tabla taxonC3mica
saveRDS(taxa, file = file.path(dir, "taxa_assignments.rds"))

library(phyloseq)
library(ggplot2)

meta <- read.csv("C:/Users/VladM/Documents/DEunir/1Cuatrim/seqyOmnics/Actividades/Actividad 3/seq_Act03WD_al03/metadatos.csv", ,
                                    row.names  = 1,
                                    fileEncoding = "UTF-16LE",
                                    stringsAsFactors = FALSE,
                                    na.strings = c("", "NA"),
                                    strip.white = TRUE) 
sample_data(ps)

# luego incorporar a phyloseq:
ps <- phyloseq(
  otu_table(seqtab.nochimP, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(meta)
)

# Transformar a abundancia relativa
ps_rela <- transform_sample_counts(ps, function(x) x / sum(x))
sample_data(ps_rela)

df_abund <- phyloseq::psmelt(ps) 
colnames(df_abund)


ggplot(df_abund, aes(x = Phylum, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Condicion.clinica, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Abundancia relativa por Phylum seg??n grupo",
       x = "Phylum",
       y = "Abundancia relativa")

ps_rela <- transform_sample_counts(ps, function(x) x / sum(x))

# Sumar por gC)nero
tax_table(ps_rela)[, "Genus"] <- as.character(tax_table(ps_rela)[, "Genus"])  # Asegura que sea texto
ps_genus <- tax_glom(ps_rela, taxrank = "Genus")  # Agrupa por Genus

# Seleccionar los 10 mC!s abundantes
top10 <- names(sort(taxa_sums(ps_genus), decreasing = TRUE))[1:8]   # de 10 a 8
ps_top10 <- prune_taxa(top10, ps_genus)
df_genus <- psmelt(ps_top10)

ggplot(df_genus, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Condicion.clinica, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 10 g??neros m??s abundantes por condici??n cl??nica",
       x = "Muestra",
       y = "Abundancia relativa")

plot_bar(ps_rela, fill = "Genus") +
  facet_wrap(~Condicion.clinica, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 50)
  ) +
  labs(
    title = "Abundancia relativa por G??nero"
  )



#Plot abundancia relativa especies

tax_table(ps)[, "Species"] <- as.character(tax_table(ps)[, "Species"])
ps_species <- tax_glom(ps_rela, taxrank = "Species")

df_all_species <- psmelt(ps_species)
top10_species <- names(sort(taxa_sums(ps_species), decreasing = TRUE))[1:8]
ps_species_top10 <- prune_taxa(top10_species, ps_species)
df_top_species <- psmelt(ps_species_top10)

ggplot(df_top_species, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Condicion.clinica, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Top 10 especies m??s abundantes por condici??n cl??nica",
       x = "Muestra",
       y = "Abundancia relativa")


# --- CURVA DE RAREFACCI??N ---

library("devtools")
library("vegan")
otu_mat <- as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)) { otu_mat <- t(otu_mat) }
rarecurve(otu_mat, step = 1000, col = "blue", cex = 0.6, label = TRUE)

# --- DIVERSIDAD ALFA Y BETA ---
# Calcular diversidad alfa.....

# Calcular todas las m??tricas de diversidad, incluyendo Simpson e InvSimpson
diversidad <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
diversidad$An??lisis <- sample_data(ps)$Condicion.clinica
# Convertir a formato largo para ggplot
library(tidyr)
library(dplyr)
diversidad_long <- pivot_longer(diversidad, cols = c("Observed", "Shannon", "Simpson", "InvSimpson"),
                                names_to = "Medida", values_to = "Valor")
# Gr??fico tipo facet con puntos
ggplot(diversidad_long, aes(x = An??lisis, y = Valor, color = An??lisis)) +
  geom_point(size = 2) +
  facet_wrap(~Medida, scales = "free_y") +
  theme_minimal() +
  labs(title = "Diversidad alfa",
       x = "An??lisis", y = "Alpha Diversity Measure") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "gray80", color = NA),
        legend.position = "right")

# --Distancia Bray-Curtis y MDS

ordu <- ordinate(ps, method = "PCoA", distance = "bray")
# Extraer porcentajes de varianza
eig_vals <- ordu$values$Relative_eig * 100
eje1 <- round(eig_vals[1], 1)
eje2 <- round(eig_vals[2], 1)

# Crear data.frame con coordenadas y metadatos
pcoa_df <- as.data.frame(ordu$vectors)
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df <- cbind(pcoa_df, as(sample_data(ps)[pcoa_df$SampleID, ], "data.frame"))

# Etiqueta personalizada (puedes usar SampleID, Condicion.clinica, etc.)
pcoa_df$Etiqueta <- pcoa_df$SampleID  # O usa: pcoa_df$Condicion.clinica

# Gr??fico con puntos y etiquetas
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Condicion.clinica, label = Etiqueta)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3.2) +  # Para etiquetas arriba de los puntos
  theme_minimal() +
  labs(title = "PCoA, distancia de Bray Curtis",
       x = paste0("Axis.1 [", eje1, "%]"),
       y = paste0("Axis.2 [", eje2, "%]"),
       color = "An??lisis") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# PERMANOVA
metadata <- data.frame(sample_data(ps))
bray_dist <- phyloseq::distance(ps, method = "bray")
adonis_result <- adonis2(bray_dist ~ Condicion.clinica, data = metadata)
print(adonis_result)


#DETECCI??N DE ASVs DIFERENCIALES ENTRE PARES DE GRUPOS
library(DESeq2)
library(phyloseq)
library(dplyr)
library(ashr)
library(DESeq2)
# Aseg??rate que Condicion.clinica sea factor
sample_data(ps)$Condicion.clinica <- factor(sample_data(ps)$Condicion.clinica)
# Convertir phyloseq a DESeq2
dds <- phyloseq_to_deseq2(ps, ~ Condicion.clinica)
# Ejecutar DESeq2
dds <- DESeq(dds, test = "Wald", fitType = "parametric")
# Comparar todos los pares posibles de condiciones cl??nicas
condiciones <- levels(sample_data(ps)$Condicion.clinica)
# Crear carpeta para guardar los resultados
res_dir <- file.path(getwd(), "deseq2_resultados")
if (!dir.exists(res_dir)) dir.create(res_dir)

for(i in 1:(length(condiciones)-1)) {
  for(j in (i+1):length(condiciones)) {
    grupo1 <- condiciones[i]
    grupo2 <- condiciones[j]
    contraste <- c("Condicion.clinica", grupo1, grupo2)
    res <- results(dds, contrast = contraste, alpha = 0.001)  # Significancia al 0.1%
    res <- lfcShrink(dds, contrast = contraste, res = res, type = "ashr")
    # Filtrar ASVs significativos
    sig_asvs <- res[which(res$padj < 0.001 & !is.na(res$padj)), ]
    sig_asvs <- sig_asvs[order(sig_asvs$padj), ]
    # Guardar resultados
    nombre_archivo <- paste0("ASVs_diferenciales_", grupo1, "vs", grupo2, ".csv")
    write.csv(as.data.frame(sig_asvs), file = file.path(res_dir, nombre_archivo))
    
    cat("Comparaci??n", grupo1, "vs", grupo2, ": ", nrow(sig_asvs), "ASVs significativos\n")
  }
}
