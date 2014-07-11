library(biomaRt)
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("hgnc_symbol")
filters <- c("chromosome_name","start","end")
genelist <- scan("excavator_results/annos/genelist.cardio.list", what = "")
myls <- list.files("excavator_results/annos/", pattern = "results.bed")

xx <- read.delim(file.path("myfile"), sep = "\t", header = F)
chr <- as.character(xx[[1]])
start <- as.character(xx[[2]])
end <- as.character(xx[[3]])
all.genes <- rep("", length(chr))
geneinlist<- rep("", length(chr))
for (j in 1: length(chr)) {
  values <- list(chromosome=chr[j],start=start[j],end=end[j])
  mygenes <- unique(getBM(attributes=attributes, filters=filters, values=values, mart=mart))
  all.genes[j] <- toString(mygenes)
}
