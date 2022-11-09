# Purpose: Analyze ASV table for beta diversity
# Date: 03/23/2022

# Load Packages ----------------------------------
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(knitr)
library(microbiome)
library(DESeq2)
library(ggplot2)
library(tools)
library(compositions)
library(plotrix)
library(psadd)
library(plyr)
library(ape)

setwd("G:/My Drive/Sow-Piglets_Elanco/R-Analysis")

## Set ggplot theme -------------------------------------------------

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}



# Load in dada2 files--------------------------------------------

load("Output/Dada2/seqtab.nochim.rda") #### RDA containing the sequence table data from dada2
load("Output/Dada2/taxa.rda") #### RDA containing the taxonomic data from dada2

# Load in sample file and build data frames--------------------------------------------
samdf <- read.csv(file = "Data/Piglets-Foglio1.csv", stringsAsFactors = TRUE)
names(samdf)[1] <- "Date"

samdf$Date <- as.character(samdf$Date)
substring(samdf$Date, 2) <- "/"

names(samdf)[2] <- "Building"
names(samdf)[3] <- "Room"
names(samdf)[4] <- "Row"
names(samdf)[5] <- "Crate"
names(samdf)[6] <- "Sex"
names(samdf)[7] <- "Sample"
samdf$Litter <- paste(samdf$Date,                             ####Concatenate into 58 litter ids
                      samdf$Building, samdf$Room,
                      samdf$Row, samdf$Crate, sep=":")

samdf$Litter <- as.factor(samdf$Litter)
samdf$LitterID <- as.numeric(factor(samdf$Litter, levels=unique(samdf$Litter)))  #### Converts concatenated litter into 58 unique litter ids


# Load in sow sample data and add SowID to samdf ---------------------------------
sows <- read.csv("Data/ElancoSowsRepr2.csv", stringsAsFactors = TRUE, header = TRUE, row.names=NULL)

sows <- subset(sows, select= -c(date1, date2, building, room, row, crate, priority, TotalBorn,
                                LiveBorn, StillBorn, Mummies, PigWean, WeanAge, W2E, PigDeaths, PWS))
samdf2 <- merge(sows, samdf, by.x = "Litter")

maltecca_metadata_piglets <- read.csv("Data/maltecca_metadata_piglets.csv", stringsAsFactors = TRUE)
names(maltecca_metadata_piglets)[1] <- "Sample"
samples.out <- rownames(seqtab.nochim)
samdf3 <- merge(maltecca_metadata_piglets, samdf2, by.x = "Sample")
rownames(samdf3) <- samples.out
samdf3$Room_Date <- as.factor(paste(samdf3$building, samdf3$room, samdf3$date2))


# Build ps object -----------------------------------------
ps_piglet <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                      sample_data(samdf3), 
                      tax_table(taxa))

# Store ASVs as short names(ASV1, ASV2, etc) -------------------------
dna <- Biostrings::DNAStringSet(taxa_names(ps_piglet))
names(dna) <- taxa_names(ps_piglet)
ps <- merge_phyloseq(ps_piglet, dna)
taxa_names(ps_piglet) <- paste0("ASV", seq(ntaxa(ps_piglet)))
ps
rm(dna, maltecca_metadata_piglets, ps_piglet, samdf, samdf2, seqtab.nochim, sows, samples.out)


## Confirm only archaea and bacteria are present-----------------------------------------------------------------------------------------------------------------------------

phyloseq::get_taxa_unique(ps, taxonomic.rank="Phylum")

## Number of taxa/species-----------------------------------------------------------------------------------------------------------------------------

phyloseq::ntaxa(ps)
## prune OTUs that are not present in any of the samples-----------------------------------------------------------------------------------------------------------------------------

ps <- phyloseq::prune_taxa(taxa_sums(ps) > 0, ps)
phyloseq::ntaxa(ps)

## Returns the total number of individuals observed from each sample-----------------------------------------------------------------------------------------------------------------------------

phyloseq::sample_sums(ps)

# Start by running some basic unconstrained ordination plots -------------
# Transform data using a VST transformation


ASV <- otu_table(ps)
min(ASV)
ASV <- ASV +1 #ensure deseq can run by adding 1 to all of data so a mean can be found
min(ASV)
ps2 <- phyloseq(otu_table(ASV, taxa_are_rows=FALSE), 
                      sample_data(samdf3), 
                      tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps2))
names(dna) <- taxa_names(ps2)
ps3 <- merge_phyloseq(ps2, dna)
taxa_names(ps2) <- paste0("ASV", seq(ntaxa(ps2)))
ps3
rm(dna,ps2)

# ps_bact <- phyloseq::subset_taxa(ps3, Kingdom=="Bacteria")
# ps_ds <- phyloseq::phyloseq_to_deseq2(ps_bact, ~1)
# ps_ds = DESeq2::estimateSizeFactors(ps_ds, type='iterate')
# ps_ds = DESeq2::estimateDispersions(ps_ds, fitType = "parametric")
# save(ps_ds, file = "Output/ps_ds.rda")

load("Output/ps_ds.rda")
png("Visualizations/plotDispEsts.png",  width = 1080, height= 670)
DESeq2::plotDispEsts(ps_ds)
dev.off()

#make a copy of the ps object with vst-transformed otu_table
ps_vst <- ps3
vst<- DESeq2::getVarianceStabilizedData(ps_ds)
phyloseq::otu_table(ps_vst) <- phyloseq::otu_table(vst, taxa_are_rows = TRUE)
ps_vst

#Now we can move into NMDS Ordination
min(otu_table(ps_vst))
ps_vst_pos <- transform_sample_counts(ps_vst, function(x) x+ 2.42)
phylum.sum = tapply(taxa_sums(ps_vst_pos), tax_table(ps_vst_pos)[, "Phylum"], sum, na.rm=TRUE)
top5phyla <- names((sort(phylum.sum, TRUE)))[1:5]
ps_vst_pos1 <- prune_taxa((tax_table(ps_vst_pos)[, "Phylum"] %in% top5phyla), ps_vst_pos)

random_tree = rtree(ntaxa(ps_vst_pos1), rooted=TRUE, tip.label=taxa_names(ps_vst_pos1))
physeq1 = merge_phyloseq(ps_vst_pos1, samdf3, random_tree)
physeq1


 # Evaluate ordination methods with Room_Date ----------------------------
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq1, dist){
  ordi = ordinate(physeq1, method=i, distance=dist)
  plot_ordination(physeq1, ordi, "samples", color="Room_Date")
}, physeq1, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=Room_Date, fill=Room_Date))
p = p + geom_point(size=4) + geom_line() 
p = p + facet_wrap(~method, scales="free")
p = p + scale_fill_brewer(type="qual", palette="Set1")
p = p + scale_colour_brewer(type="qual", palette="Set1")
png("Visualizations/ord_Rm_dt.png",  width = 1080, height= 670)
p
dev.off()

# Evaluate ordination methods with Litter ----------------------------
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist2 = llply(as.list(ord_meths), function(i, physeq1, dist){
  ordi = ordinate(physeq1, method=i, distance=dist)
  plot_ordination(physeq1, ordi, "samples", color="Litter")
}, physeq1, dist)
names(plist2) <- ord_meths
pdataframe2 = ldply(plist2, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe2)[1] = "method"

p2 = ggplot(pdataframe2, aes(Axis_1, Axis_2, color=Litter, fill=Litter))
p2 = p2 + geom_point(size=4) + geom_line() 
p2 = p2 + facet_wrap(~method, scales="free")
p2 = p2 + scale_fill_brewer(type="qual", palette="Set1")
p2 = p2 + scale_colour_brewer(type="qual", palette="Set1")
png("Visualizations/ord_litter.png",  width = 1080, height= 670)
p2
dev.off()



# ord1 <- phyloseq::ordinate(physeq1, "MDS", "bray", weighted= TRUE, trymax=100)
# # Evaluate NMDS run with stress
# 
# 
# p1 <- phyloseq::plot_ordination(ps_vst_pos1, ord1,
#                           type="taxa",
#                           color="Phylum",
#                           title="MDS Bray Taxa") +
#   ggplot2::geom_line()  +
#   ggplot2::geom_point(size=5) +
#   facet_wrap (~Phylum, 3)
# p1

# Heat map of Litter
(p <- plot_heatmap(ps, "Pearson", "bray", "Litter", "Family", 
                   low= "#FFFFCC", high= "#000033", na.value="white"))
heatmap(otu_table(ps_vst_pos1))
