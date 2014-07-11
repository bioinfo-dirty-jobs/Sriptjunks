library("GenomeGraphs")
library("biomaRt")
setwd("/home/maurizio/Desktop/")

nome = "DDIT3"

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# fetch HGNC symbols chromosome for 
hgnc <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id","chromosome_name"), filters = "hgnc_symbol",value=nome, mart = mart)

gene <-makeGene("hgnc_symbol", id = nome , biomart = mart)
gene.transcripts <- makeTranscript("hgnc_symbol", id = nome , biomart = mart)
# make an axis 5'-3' for (+), 3'-5' for (-)
gene.strand <- as.character(gene@ens$strand[1])
if(gene.strand ==  "1") {chrom.axis <- makeGenomeAxis(add53 = TRUE)} else
  if(gene.strand == "-1") {chrom.axis <- makeGenomeAxis(add35 = TRUE)}
poli <- list(makeTitle(nome), makeIdeogram(chromosome = hgnc$chromosome_name[1]),  chrom.axis, "gene" = gene, "transcripts" = gene.transcripts)
####################################################################################################################à
png(filename = "Prova_out.png", width = 800, height = 600)
gdPlot(poli)
dev.off()

exit()
