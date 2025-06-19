
#DADA2
###Instalar librerias

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("dada2", quietly = TRUE))
  BiocManager::install("dada2", version = "3.21")
####Cargar librerias

library("dada2")
library("stats")
library("BiocParallel")
library("R.utils")
####DADA2 PIPELINE 2.0 - DATASET PROYECTO FUSARIUM- Revisado Octubre 2024

#Establecer ruta a archivos limpios
setwd("C:/Users/PC/Documents/INIAP/Metagenomica/Tejido/FASTQ_Tejido/Bacterias/Tejido_16S_Fastq_Comprimida")
path <- "C:/Users/PC/Documents/INIAP/Metagenomica/Tejido/FASTQ_Tejido/Bacterias/Tejido_16S_Fastq_Comprimida"
list.files(path)

#Los archivos forward y reverse tienen el formato: _R1_001.fastq y _R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) #para correr despues de trimmomatic
print(fnFs)
print(length(fnFs))

fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
print(fnRs)
R.utils::countLines(fnFs[1])

#Extraer los nombres de las muestras 
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
print(sample.names)

print(fnFs)
print(length(fnFs))

#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
file.exists(fnFs[1])
dev.new()
plotQualityProfile(fnFs[1])
plotQualityProfile(fnFs[1:2])


# Trimming
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names

names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)                  #for ITS sequencing consider leave out truncLen
#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Plot error
plotErrors(errF, nominalQ=TRUE)

##Dereplicate (this step does not be included in version 1.16 but it is in version 1.8)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
save.image(file="/ArchivosMC/Documentos/ITS_/Día1hastaderep.rda")

##Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)  #continuando el proceso luego de la dereplicacion añadir pool=TRUE
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

save.image(file="/ArchivosMC/Documentos/ITS_M/Día1hastasampleinference.rda")


dadaPFs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)  
dadaPRs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

#Inspect returned class object
dadaFs[[9]]
dadaRs[[9]]
dadaPFs[[9]]
dadaPRs[[9]]

#Merge paired reads to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)  #con dereplicacion
mergersP <- mergePairs(dadaPFs, derepFs, dadaPRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[9]])
head(mergersP[[9]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
seqtabP <- makeSequenceTable(mergersP)

## The sequences being tabled vary in length.
dim(seqtab)
dim(seqtabP)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
table(nchar(getSequences(seqtabP)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

seqtab.nochimP <- removeBimeraDenovo(seqtabP, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochimP)
sum(seqtab.nochimP)/sum(seqtabP)

#Track reads trough pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

getN <- function(x) sum(getUniques(x))
trackP <- cbind(out, sapply(dadaPFs, getN), sapply(dadaPRs, getN), sapply(mergersP, getN), rowSums(seqtab.nochimP))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track

colnames(trackP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(trackP) <- sample.names
head(trackP)
trackP

#saveRDS(track, "/ArchivosMC/Documentos/ITS_MENDOZA2/summary_table18S.rds")
saveRDS(track, "/ArchivosMC/Documentos/ITS_MENDOZA2/summary_table18S.rds")
save.image(file="/ArchivosMC/Documentos/ITS_MENDOZA2/summary_table18S.rda")

###ASSIGN TAXONOMY### (corrió la de la carpeta databases de todos)
taxa <- assignTaxonomy(seqtab.nochim, "/ArchivosMC/Descargas/UNITEFUNGI.fasta", multithread=TRUE)

#ASSIGN SPECIES STEP
taxa <- addSpecies(taxa, "/home/esteban.mena/databases/silva_nr_v138_train_set.fa.gz", verbose=TRUE)

#Inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

save.image(file="/home/taxonomicassigments")
saveRDS(seqtab.nochim, "/home/18S_M/seqtab.nochim2.rds")
saveRDS(taxa, "/home/18S_M/taxa2.rds")
save.image(file="/home/18S_ME/dada22.rda")

#Phyloseq
#Instalación de paquetes comparación
install.packages('pacman')
BiocManager::install("phyloseq")
library(pacman)
packages <- c('phyloseq', 'Biostrings', 'ggplot2', 'tidyverse', 'vegan', 'decontam', 'ape', 'DESeq2', 'microbiome', 'metagenomeSeq', 'remotes')
pacman::p_load(char = packages)
remotes::install_github("MadsAlbertsen/ampvis2", force=TRUE)
remotes::install_github("adw96/breakaway", force=TRUE)
BiocManager::install("limma")
BiocManager::install("DECIPHER")
install.packages("phangorn")
BiocManager::install("tidyverse")
library(tidyverse)
library(limma)
library(phyloseq); packageVersion("phyloseq")
library(phangorn)
library(dada2)
library(DECIPHER)
library(ggplot2); packageVersion("ggplot2")
library(parallel)
#Temas (colores) de ggplot
theme_set(theme_bw())

samples.out <- rownames(seqtab.nochim)
samples.out
sitio <- sapply(strsplit(samples.out, "M"), `[`, 1)

tratamiento <- substr(sitio,1,1)
#sitio <- substr(sitio,2,999)
muestra <- sapply(strsplit(samples.out, "M"), `[`, 2)
tejido <- substr(muestra,1,1)
replica <- substr(muestra,2,2)
samdf <- data.frame(Tratamiento=tratamiento, Tejido=tejido, Réplica=replica)
samdf$Análisis <- "Deflagración"
samdf$Análisis[samdf$Tratamiento<2] <- "Remanente"
samdf$Análisis[samdf$Tratamiento>2] <- "Predetonación"
rownames(samdf) <- samples.out
table

#Almacenar tabla de metadatos en outdir
write.csv(samdf, file = file.path(path, 'metadatos2.csv'))

save.image(file="/ArchivosMC/Descargas/ITS/tablametadatos.rda")


#Instalación de paquetes
packages <- c('phyloseq', 'Biostrings', 'ggplot2', 'tidyverse', 'vegan', 'decontam', 'ape', 'DESeq2', 'microbiome', 'metagenomeSeq', 'remotes')

pacman::p_load(char = packages)
remotes::install_github("MadsAlbertsen/ampvis2", force=TRUE)

install.packages("remotes")
remotes::install_github("kasperskytte/ampvis2", Ncpus = 6, force=TRUE)

remotes::install_github("adw96/breakaway")
BiocManager::install("limma", force=TRUE)
BiocManager::install("DECIPHER")
install.packages("phangorn")
BiocManager::install("tidyverse")
library(tidyverse)
library(limma)
library(phyloseq)
library(phangorn)
library(dada2)
library(DECIPHER)

#Crear un archivo fasta con los ASVs
#Preparar secuencias y ancabezados
secuencias <- colnames(seqtab.nochim)
encabezados <- paste('>ASV', 1:ncol(seqtab.nochim), sep = '_')
#Interlinear
asv_fasta <- c(rbind(encabezados, secuencias))

#Preparar secuencias
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs 

#Alinear
alignment <- AlignSeqs(DNAStringSet(seqs),
                       anchor = NA,
                       iterations = 5,
                       refinements = 5)
phang.align <- phyDat(as(alignment, 'matrix'), type = 'DNA')
dm <- dist.ml(phang.align)

treeNJ <- NJ(dm) 
fit = pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR,
                    model = 'GTR',
                    optInv = TRUE,
                    optGamma = TRUE,
                    rearrangement = 'stochastic',
                    control = pml.control(trace = 0))

#alternativamente guardar objeto fitGTR que contiene el arbol y otros objetos generados
saveRDS(fitGTR, '/home/18S_MENDOZA/fitGTR2.rds')
saveRDS(alignment, '/home/18S_MENDOZA/alignment2.rds')
saveRDS(phang.align, '/home/18S_MENDOZA/phang.align2.rds')
saveRDS(treeNJ, '/home/18S_MENDOZA/treeNJ2.rds')
saveRDS(fit, '/home/18S_MENDOZA/fit2.rds')
load(file="/ArchivosMC/Descargas/taxonomyassignabundance.rda")
#elemento grande 
saveRDS(dm, '/home/dm2.rds')

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa), phy_tree(fitGTR$tree))
#Almacenar el objeto de phyloseq en el outdir
saveRDS(ps, file = file.path(outdir, 'ps_18S2.rds'))
saveRDS(ps, '/home/18S_MENDOZA/ps2_18S.rds')

#Cambiar el nombre de los ASVs en el objeto ps
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

# Juntar el objeto de secuencias al objeto de phyloseq:
ps<- merge_phyloseq(ps, dna)
#Renombrar ASVs
taxa_names(ps) <- paste("ASV", 1:ntaxa(ps), sep = "_")
#Observar el cambio en los primeros 3 ASVs
taxa_names(ps)[1:3]

#Eliminar taxones eucariotas
# Crear subsets para cloroplast, mitocondria, and eukaryotes (solo quité las bacterias, porque cloroplasto y mitocondria no hay):
chlr <- subset_taxa(ps, Order == "Cloroplasto")
mit <- subset_taxa(ps, Family == "Mitocondria")

bac <- subset_taxa(ps, Kingdom == "Bacteria")

# Generar una variablew con todos los taxones a remover:
elim_tax <- c(taxa_names(bac))
# Mantener taxones deseados:
taxa <- taxa_names(ps)
des_tax <- taxa[!(taxa %in% elim_tax)]

# Objeto ps reducido conteniendo solo subset deseado:
ps_red <- prune_taxa(des_tax, ps)

#Determinar proporción que se mantuvo en ps_red
pre <- sum(sample_sums(ps))
post <- sum(sample_sums(ps_red))
post / pre

#Remover muestras con menos de 1000 conteos (REVISAR ERROR INVALID CLASS SAMPLE DATA MUST HAVE ZERO DIMENSIONS)
ps2 <- subset_samples(ps_red, sample_sums(ps_red) > 1000)
sample_sums
library(phyloseq)



#Observar datos en el objeto
sample_data(ps)[1:9]
head(t(otu_table(ps)[1:9]))
head(tax_table(ps2))

update.packages("phyloseq")

#Normalizacion
# Normalizacion por proporcion:
library(microbiome)
ps_prop <- transform(ps, "compositional")
# Observar tabla de OTUs:
otu_table(ps_prop)[1:9, 1:9]
# Observar el conteo por muestra, todos deben sumar 1:
head(sample_sums(ps_prop))

#Graficar abundancia sin normalizar
# Tomar subset de Remanente (Análisis Remanente)
ps_prop <- subset_samples(ps_red, Análisis = "Remanente")
ps_prop1 <- subset_samples(ps_red, Análisis = "Deflagración")
ps_prop2 <- subset_samples(ps_red, Análisis = "Predetonación")

# Agglomerar a nivel de Phylum:
pseud_phylum <- tax_glom(ps_prop, taxrank = "Phylum", NArm = FALSE)
pseud_phylum1 <- tax_glom(ps_prop1, taxrank = "Phylum", NArm = FALSE)
pseud_phylum2 <- tax_glom(ps_prop2, taxrank = "Phylum", NArm = FALSE)
# Remover phyla poco comun:
pseud_phylum <- subset_taxa(pseud_phylum, taxa_sums(pseud_phylum) > 0.1)
pseud_phylum1 <- subset_taxa(pseud_phylum1, taxa_sums(pseud_phylum1) > 0.1)
pseud_phylum2 <- subset_taxa(pseud_phylum2, taxa_sums(pseud_phylum2) > 0.1)
# Graficar:
plot_bar(pseud_phylum, fill = "Phylum") +
  facet_wrap(~Análisis, scales = "free_x", nrow = 1)

save.image(file="/ArchivosMC/Descargas/hastafitGTR.rda")

load(file="/ArchivosMC/Descargas/mary/ITS_ME/hastaabundanciataxonómicanormalizadaphylumygenero.rda")
