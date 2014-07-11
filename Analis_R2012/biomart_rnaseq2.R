library(biomaRt)
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("hgnc_symbol","chromosome_name","chromosome_location")
#filters <- c("chromosome_name","start","end")
filters=c("refseq_mrna")

#mygenes <- unique(getBM(attributes=attributes, filters= filters, values=c("NM_001008701"),mart=mart))
#grep("start",attributi, ignore.case=T, value=T)