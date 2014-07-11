ibrary(biomaRt)
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


#--------------------------------------
library("AnnotationForge")
ls("package:GO.db")
ls("package:AnnotationForge")
library(GOstats)
# Specify the parameters for the enrichment analysis
params = new("GOHyperGParams", # Type of object to              be created
             geneIds = de.Entrez.ids, # gene set on which             enrichment will be performed
             annotation  = "org.Hs.eg.db", #   database with   annotation
             ontology = "BP", # GO to look for
             pvalueCutoff = 0.05, # p-value threshold
             conditional = FALSE, #  a logical indicating    whether the calculation should condition on the GO structure 
             testDirection ="over")
# looking for over-  represented terms
# Perform the hypeGeometric test with the above-defined parameters
result = hyperGTest(params)
# Store the results into the 'Go_enrichment.html'file under the 'deseq' folder
htmlReport(result, "resulst/GO_enrichment.html")

             