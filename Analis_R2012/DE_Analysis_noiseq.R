# DE Analysis with NOISeq R package
## load libraries and data
library(NOISeq) 
library(goseq)
library(GO.db)
data(Marioni)
head(mycounts) 

myfactors = data.frame(Tissue = c("Kidney", "Liver", "Kidney", "Liver", "Liver", "Kidney", "Liver", "Kidney", "Liver", "Kidney"), TissueRun = c("Kidney_1", "Liver_1", "Kidney_1", "Liver_1", "Liver_1", "Kidney_1", "Liver_1", "Kidney_2", "Liver_2", "Kidney_2")) 
head(mylength)
head(mybiotypes) 

head(mychroms)

## import table of counts in R  

mydata <- readData(data = mycounts, length = mylength, biotype = mybiotypes, chromosome = mychroms, factors = myfactors)
mydata <- readData(data = mycounts, chromosome = mychroms, factors = myfactors) 
mydata <- addData(mydata, length = mylength, biotype = mybiotypes) 

## Generating data for exploratory plots
myexplodata <- dat(mydata, type = "biodetection", selection = c(3, 4)) 
## ## Biodetection
mybiodetection <- dat(mydata, selection = c(1, 2), k = 0, type = "biodetection")
explo.plot(mybiodetection) 
## ## Count distribution comparison
mycd <- dat(mydata, selection = c(1, 2), type = "cd")
explo.plot(mycd) 
## ## Distribution of counts per biological group
mycountsbio <- dat(mydata,selection=NULL,k=0,ndepth=5,type="countsbio") 
explo.plot(mycountsbio,toplot=1,samples=c(1,2),ylim=c(0,200)) 
explo.plot(mycountsbio,toplot="protein_coding",samples=5,ylim=c(0, 200)) 
## ## Length of detected features
myDLbio <- dat(mydata, selection = NULL, k = 0, ndepth = 5, type = "DLbio") 
explo.plot(myDLbio, samples = 3:6, toplot = "protein_coding") 
## ## Saturation plots 
mysaturation <- dat(mydata, k = 0, ndepth = 5, newdetections = TRUE, type = "saturation") 
explo.plot(mysaturation, toplot = 1, samples = 1:2, ylim = NULL, yrightlim = NULL) 
explo.plot(mysaturation, toplot = "protein_coding", samples = 1:4, ylim = NULL, yrightlim = NULL) 


## Differential expression
## ## normalization 
myRPKM=rpkm(assayData(mydata)$exprs,long=mylength,k=0,lc=1) 
myUQUA=uqua(assayData(mydata)$exprs,long=mylength,lc=0.5,k=0)
myTMM=tmm(assayData(mydata)$exprs,long=1000,lc=1) 
head(myRPKM[,1:4]) 
## ## DE testing (NOISeq-real: using available replicates)
mynoiseq <- noiseq(mydata, k = 0.5, norm = "rpkm", factor = "Tissue", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "technical") 
head(mynoiseq@results[[1]])
mynoiseq.tmm = noiseq(mydata, k = 0.5, norm = "tmm", factor = "TissueRun", conditions = c("Kidney_1", "Liver_1"), lc = 0, replicates = "technical") 
## ## DE testing (NOISeq-sim: no replicates available ) 
myresults <- noiseq(mydata, factor = "Tissue", k = NULL, norm = "n", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no") 

## Output data
mytable <- mynoiseq.tmm@results[[1]]
head(mytable)

## How to select the differentially expressed features
mynoiseq.deg <- degenes(mynoiseq, q = 0.8, M = NULL) 
head(mynoiseq.deg)
mynoiseq.deg1 <- degenes(mynoiseq, q = 0.8, M = "up")
mynoiseq.deg2 <- degenes(mynoiseq, q = 0.8, M = "down")
DE.plot(mynoiseq, q = 0.8, graphic = "MD", ylim = c(0, 50)) 
DE.plot(mynoiseq, chromosomes = c("1", "2"), log.scale = TRUE, join = FALSE, q = 0.8, graphic = "chrom") 
##pdf("manhattan.pdf", width = 12, height = 50) 
##DE.plot(mynoiseq, chromosomes = c("I","II"), log.scale = TRUE, join = FALSE, q = 0.8) 
##dev.off() 


genes.prob <- mytable$prob
genes.id <- rownames(mytable)
genes.prob[which(genes.prob >= 0.8)] <- 1
genes.prob[which(genes.prob < 0.8)] <- 0
genes.prob[which(is.na(genes.prob))] <- 0

names(genes.prob) <- genes.id

## write the genelist of DE genes
cat(genes.id[which(genes.prob == 1)], file = "~/mycalls/genelist.txt", sep = "\n")
##
## perform GSE analysis with GOseq
par(mfrow = c(1,1))
pwf <- nullp(genes.prob,'hg19','ensGene')
GO.pvals <- goseq(pwf,'hg19','ensGene')
GO.pvals.filtered <- GO.pvals[which(GO.pvals$over_represented_pvalue < 0.001), ]
xx <- as.list(GOTERM)
ncalls <- dim(GO.pvals.filtered)[1]
myterms <- rep("", ncalls)
mydefinitions <- rep("", ncalls)
myontologies <- rep("", ncalls)
for (i in 1: ncalls) {
  myterms[i] <- Term(xx[[GO.pvals.filtered[i,1]]])
  mydefinitions[i] <- Definition(xx[[GO.pvals.filtered[i,1]]])
  myontologies[i] <- Ontology(xx[[GO.pvals.filtered[i,1]]])
}
GO.calls <- cbind(GO.pvals.filtered, myterms, mydefinitions, myontologies)
colnames(GO.calls) <- c(names(GO.pvals.filtered), "term", "def", "Onto")
write.table(GO.calls, file = "~/mycalls/GO.results.txt.csv", sep = "\t", row.names = F, quote = F)








