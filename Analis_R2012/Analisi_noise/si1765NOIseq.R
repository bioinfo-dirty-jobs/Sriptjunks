#analisi 1765
library(NOISeq) 
library(goseq)
library(GO.db)
#set folder path
setwd('/home/maurizio/Analis_R2012/Analisi_noise/')
a<-read.table("Dati_countsi1765.csv",header=T,sep="\t")
names(a)<-c("row.names","siGlo","siGlo","siGlo","siFC","siFC","siFC")
b<-a[,-1]
row.names(b)<-a[,1]
b[is.na(b)] <- 0 ## sostituisci NA con 0
#carica dati
myfactor=data.frame(Corsa=c("siGlo","siGlo","siGlo","siFC","siFC","siFC"),monica=c("A","A","A","B","B","B"))
mydata <- readData(data = b,  factors = myfactor)


                 
mynoiseq <- noiseq(mydata, k = 0.5, norm = "rpkm", factor = "Corsa", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "technical") 
mynoiseq.tmm = noiseq(mydata, k = 0.5, norm = "tmm", factor = "Corsa", conditions = c("siGlo", "siFC"), lc = 0, replicates = "biological") 
## ## DE testing (NOISeq-sim: no replicates available ) 
png("Plot_tmm.png")
DE.plot(mynoiseq.tmm, q = 0.8, graphic = "MD")
dev.off()
png("Plot_rpkm.png")
DE.plot(mynoiseq, q = 0.8, graphic = "MD", ylim = c(0, 50))
dev.off()


pdf("manhattan.pdf", width = 12, height = 50)
DE.plot(mynoiseq.tmm, chromosomes = c("I","II"), log.scale = TRUE, join = FALSE, q = 0.8)
dev.off()

head(mynoiseq)
mytable_tm <- mynoiseq.tmm@results[[1]]
head(mytable)
mytable_fpkm <- mynoiseq@results[[1]]
head(mytable)
write.table(mytable_tm,"noiseqresults_tmm.csv",sep="\t")
write.table(mytable_fpkm,"noiseqresults_fpkm.csv",sep="\t")
## How to select the differentially expressed features
mynoiseq.deg <- degenes(mynoiseq, q = 0.8, M = NULL) 
mynoiseq.deg_tmm <- degenes(mynoiseq.tmm, q = 0.8, M = NULL) 
mynoiseq.deg1 <- degenes(mynoiseq, q = 0.8, M = "up")
mynoiseq.deg2 <- degenes(mynoiseq, q = 0.8, M = "down")

mynoiseq.deg1_tmm <- degenes(mynoiseq.tmm, q = 0.8, M = "up")
mynoiseq.deg2_tmm <- degenes(mynoiseq.tmm, q = 0.8, M = "down")
#writeUp

write.table(mynoiseq.deg1_tmm,"Geni_up_tmm.csv",sep="\t")
write.table(mynoiseq.deg1,"Geni_up_fpkmm.csv",sep="\t")

