library(Biobase)
library(GEOquery)
library(limma)
library("chimera")
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
fusion.names <- fusionName(tmp)
fusion.names
myset <- tmp[[13]]
tmp.seq <- chimeraSeqs(myset, type="transcripts")
#write.XStringSet(tmp.seq, paste(sub(":","_",fusion.names[[13]]),".fa",sep=""), format="fasta")
library("chimera")
setwd("/home/maurizio/Desktop/GeniFusioneLeio/Pazienti_defuse/")

tmp <- importFusionData(format="defuse","A1a__def_results.filtered.tsv")
fusion.names <- fusionName(tmp)
fusion.names