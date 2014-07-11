#proa noiseq
setwd('/home/maurizio/Desktop/NGS_work/NOISeqBIO_v00/CountExperiments_Rnaseq/')
system("ls")
A <-read.table("SampleA.count.txt")
B <-read.table("SampleB.count.txt")
C <-read.table("SampleC.count.txt")
D <-read.table("SampleD.count.txt")
E<-read.table("SampleE.count.txt")
fare <-read.table("SampleF.count.txt")
countA<-as.numeric(as.character(A[,1]))
countB<-as.numeric(as.character(B[,1]))
countC<-as.numeric(as.character(C[,1]))
countD<-as.numeric(as.character(D[,1]))
countE<-as.numeric(as.character(E[,1]))
countFA<-as.numeric(as.character(Fa[,1]))

lista1=merge(A,B,by="V2",all=T)
lista2=merge(lista1,C,by="V2",all=T)
lista3=merge(D,E,by="V2",all=T)
lisata4 =merge(lista3,fare,by="V2",all=T)
Tutto<-merge(lista2,lisata4,by="V2",all=T)
names(Tutto)<-c("gene_names",siGlo","siGlo","siGlo","siFC","siFC","siFC")
write.table(Tutto,"Dati_countsi1765.csv",sep="\t",row.names=F)
#analisi 1765
library(NOISeq) 

a<-read.table("Dati_countsi1765.csv",header=T,sep="\t",row.names=)
names(a)<-c("row.names","siGlo","siGlo","siGlo","siFC","siFC","siFC")

myfactors = data.frame(Tissue = c("siGlo", "siFC"), Replicate = c("SiGlo_1", "SiGlo_2", "SiGlo_3", "siFC_1", "siFC_2", "siFC_3")) 

mydata <- readData(data = a[,2:6],  factors = myfactors)

names(b)<-c("siGlo","siGlo","siGlo","siFC","siFC","siFC")
mydata<-readData(data=b,factors=myfactors)
factors = data.frame( Situ=c("siGlo1","siGlo2","siGlo3","siFC1","siFC2","siFC3")

b[is.na(b)] <- 0



mynoiseq <- noiseq(mydata, k = 0.5, norm = "rpkm", factor = "Cond", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "technical") 
mynoiseq.tmm = noiseq(mydata, k = 0.5, norm = "tmm", factor = "Cond", conditions = c("siGlo", "siFC"), lc = 0, replicates = "biological") 
## ## DE testing (NOISeq-sim: no replicates available ) 
DE.plot(mynoiseq.tmm, q = 0.8, graphic = "MD", ylim = c(0, 50))
pdf("DEplotmm.pdf")
DE.plot(mynoiseq.tmm, q = 0.8, graphic = "MD", ylim = c(0, 50))
dev.off()
pdf("manhattan.pdf", width = 12, height = 50)
DE.plot(mynoiseq.tmm, chromosomes = c("I","II"), log.scale = TRUE, join = FALSE, q = 0.8)
dev.off()
mynoiseq <- noiseq(mydata, k = 0.5, norm = "rpkm", factor = "Cond", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "technical")
head(mynoiseq)
mytable <- mynoiseq.tmm@results[[1]]
head(mytable)
mytable <- mynoiseq@results[[1]]
head(mytable)
mynoiseq.tmm = noiseq(mydata, k = 0.5, norm = "tmm", factor = "Cond", conditions = c("siGlo", "siFC"), lc = 0, replicates = "biological")
mytable2 <- mynoiseq.tmm@results[[1]]
head(mytable2)
write.table(mytable2,"noiseqresults.csv",sep="\t")
write.table(mytable2,"noiseqresults_tmm.csv",sep="\t")
write.table(mytable,"noiseqresults_rpkm.csv",sep="\t")
sytem("rm noiseqresults.csv")
sytems("rm noiseqresults.csv")
system("rm noiseqresults.csv")

