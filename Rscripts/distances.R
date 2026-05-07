library(Rcpp)
library(dada2)
install.packages("BiocManager")
BiocManager::install("dada2")
#secuencias-calidad
path <- "C:/Users/juanj/OneDrive/Escritorio/Abejas-Microbiota/Abejas-Microbiota/Analisis"
list.files(path)
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`,1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Filtrar
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
#maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

#Errores
errF <- learnErrors(filtFs, multithread = T)
errR <- learnErrors(filtRs, multithread = T)
plotErrors(errF, nominalQ=TRUE)

#sample inference 
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)
dadaFs[[1]]

#merge 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

#sequence table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)\sum(seqtab)
write.csv(seqtab.nochim,"/datacnmat01/biologia/gimur/proyectos/microbiomas/7_25_Abejas/Proyecto_Andre/Data_Analysis/Analisis/seqtab_nochim_bee.csv")


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), 
               #getN(dadaFs), getN(dadaRs), getN(mergers), 
               rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "/datacnmat01/biologia/gimur/proyectos/microbiomas/7_25_Abejas/Proyecto_Andre/Data_Analysis/Analisis/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/datacnmat01/biologia/gimur/proyectos/microbiomas/7_25_Abejas/Proyecto_Andre/Data_Analysis/Analisis/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print
tail (taxa.print)
write.csv(taxa,"/datacnmat01/biologia/gimur/proyectos/microbiomas/7_25_Abejas/Proyecto_Andre/Data_Analysis/Analisis/taxa_bee_16s.csv")
#transponer tabla sin quimeras 

p <- t(read.csv("/datacnmat01/biologia/gimur/proyectos/microbiomas/7_25_Abejas/Proyecto_Andre/Data_Analysis/Analisis/seqtab_nochim_bee.csv"))
write.csv(p,"/datacnmat01/biologia/gimur/proyectos/microbiomas/7_25_Abejas/Proyecto_Andre/Data_Analysis/Analisis/seqtab_nochim_transp.csv")


