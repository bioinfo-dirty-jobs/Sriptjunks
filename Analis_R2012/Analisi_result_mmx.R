library(biomaRt)
setwd('~/Desktop/B_analisi/')
data<-read.table("B_result.mmseq",head=T,sep="\t",stringsAsFactors=F)
db2query = useMart("ensembl")
listDatasets(db2query)
#grep("homo",listDatasets(db2query))
mart <- useDataset("hsapiens_gene_ensembl",db2query)
#pages = attributePages(ensembl)
#istAttributes(mart,page="feature_page") search information on typem of list files


ensembl2entrez = getBM(filters="ensembl_transcript_id", attributes =c("ensembl_gene_id", "entrezgene","unigene"), values =data[,1], mart = mart)
de.Entrez.ids =  unique(ensembl2entrez$entrezgene[!is.na(ensembl2entrez$entrezgene) ] )
library("GO.db")

biocLite("GO.db")


library("AnnotationForge")
ls("package:GO.db")
ls("package:AnnotationForge")
library(GOstats)
library("org.Hs.eg.db")
# Specify the parameters for the enrichment analysis
params = new("GOHyperGParams", geneIds = de.Entrez.ids, # gene set on which             enrichment will be performed
             annotation  = "org.Hs.eg.db", #   database with   annotation
             ontology = "BP", # GO to look for
             pvalueCutoff = 0.05, # p-value threshold
             conditional = FALSE, #  a logical indicating    whether the calculation should condition on the GO structure 
             testDirection ="over")
# looking for over-  represented terms
# Perform the hypeGeometric test with the above-defined parameters
result = hyperGTest(params)
# Store the results into the 'Go_enrichment.html'file under the 'deseq' folder
htmlReport(result, "results/GO_enrichment.html")
#prova2
params2 = new("GOHyperGParams", geneIds = de.Entrez.ids, # gene set on which             enrichment will be performed
             annotation  = "org.Hs.eg.db", #   database with   annotation
             ontology = "BP", # GO to look for
             pvalueCutoff = 0.01, # p-value threshold
             conditional = FALSE, #  a logical indicating    whether the calculation should condition on the GO structure 
             testDirection ="over")
# looking for over-  represented terms
# Perform the hypeGeometric test with the above-defined parameters
result2 = hyperGTest(params2)
# Store the results into the 'Go_enrichment.html'file under the 'deseq' folder
htmlReport(result2, "results/GO_enrichment2.html")


