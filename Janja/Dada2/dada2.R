## ----libraries----------------------------
library(dada2)

##-----load in config file ---------------
source("Dada2_Config.R")

out_dir <- "Dada2_output"
Janja_master_dir <- Janja_master_dir
Multiplexed_Seqs_Directory <- Multiplexed_Seqs_Directory
path <- paste0(sub("/[^/]+$", "", Multiplexed_Seqs_Directory),"/Demultiplexed_Seqs")
set_maxN <- maxN
set_maxEE <- maxEE
f_fastq_format <- "_R1_001.fastq"
r_fastq_format <- "_R2_001.fastq"
wd <- getwd()

# ---Prepare directories ------------------------------
ifelse(!dir.exists(file.path(wd, out_dir)), dir.create(file.path(wd, out_dir)), FALSE)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern=f_fastq_format, full.names = TRUE))
fnRs <- sort(list.files(path, pattern=r_fastq_format, full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## ----qualityplot------------------------------------------
fnFs_path <- file.path(out_dir,"fnFs_quality.png")
png(filename= fnFs_path, width=1000, height=500)
plotQualityProfile(fnFs[1:2])
dev.off()

## ----quality plot-----------------------------------------
fnRs_path <- file.path(out_dir,"fnRs_quality.png")
png(filename = fnRs_path, width=1000, height=500)
plotQualityProfile(fnRs[1:2])
dev.off ()

## Path to filtered fastq files---------------------------------------------------------
filtFs <- file.path(paste0(sub("/[^/]+$", "", Multiplexed_Seqs_Directory)), "Filtered_Seqs", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(paste0(sub("/[^/]+$", "", Multiplexed_Seqs_Directory)), "Filtered_Seqs", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


## ----filter and trim raw fastq files--------------------------------------
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              truncLen=c(f_truncLen,r_truncLen), maxEE=set_maxEE, 
              maxN=set_maxN, rm.phix=TRUE, trimLeft = c(f_primer_len, r_primer_len),
              compress=TRUE, verbose=TRUE, matchIDs=TRUE, multithread=TRUE)

## write the output for tracking---------------------------------------------------------
out_df <- as.data.frame(out)
out_df$Retained <- out_df$reads.out/out_df$reads.in
head(out_df)

## ----learn error------------------------------------------
errF <- learnErrors(filtFs, multithread=TRUE)
errF_path <- file.path(out_dir,"errF.rda")
save(errF, file= errF_path)

## ----learn error2-----------------------------------------
errR <- learnErrors(filtRs, multithread=TRUE)
errR_path <- file.path( out_dir,"errR.rda")

save(errR, file = errR_path)

## ----plot error-------------------------------------------
plot_errF_path <- file.path(out_dir,"plot_errF.png")
png(filename= plot_errF_path, width=1000, height=500)
plotErrors(errF, nominalQ=TRUE)
dev.off()

## ----plot error2------------------------------------------
plot_errR_path <- file.path(out_dir,"plot_errR.png")
png(filename= plot_errR_path, width=1000, height=500)
plotErrors(errR, nominalQ=TRUE)
dev.off()

## ---- dereplicate fastq files--------------------------------------
derepFs <- derepFastq(filtFs, verbose=TRUE)

derepRs <- derepFastq(filtRs, verbose=TRUE)
#  Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


## View first sample in dadaFs---------------------------------------------------------
dadaFs[[1]]


## ----merger, results='hide'-------------------------------
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

## ---------------------------------------------------------
head(mergers[[1]])

## ----seqtab-----------------------------------------------
seqtab <- makeSequenceTable(mergers)

## ---------------------------------------------------------
dim(seqtab)


## ---------------------------------------------------------
table(nchar(getSequences(seqtab)))


## ---------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
dim(seqtab.nochim)

seqtab.nochim_path <- file.path(sub("/[^/]+$", "", Multiplexed_Seqs_Directory),"ASV_table.csv")

write.csv(seqtab.nochim, file=seqtab.nochim_path)

## ---------------------------------------------------------
sum(seqtab.nochim)/sum(seqtab)


## ---------------------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", 
                     "denoisedR", "merged", "nonchim")

track <- as.data.frame(track)

track$Tot_Perc_Retained <- track$nonchim/track$input

head(track)

track_path <- file.path(out_dir,"track_seq_retention.csv")
write.csv(track, file= track_path)

## ----assign taxa------------------------------------------
train_destfile= paste0(wd,"/silva_nr99_v138.1_train_set.fa.gz")
train_fileURL <- "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1"
 if (!file.exists(train_destfile)) {
    download.file(train_fileURL ,train_destfile,method="auto") }
 
species_destfile= paste0(wd,"/silva_species_assignment_v138.1.fa.gz")
species_fileURL <- "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1"
 if (!file.exists(species_destfile)) {
    download.file(species_fileURL ,species_destfile,method="auto") }

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", 
                       multithread=TRUE)

taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

## ---------------------------------------------------------
taxa_path <- file.path(sub("/[^/]+$", "", Multiplexed_Seqs_Directory),"taxa.csv")

write.csv(taxa, file=taxa_path)

## ---------------------------------------------------------
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

print(paste("--------------Dada2 Pipeline Finished", Sys.time(),"-----------------------"))