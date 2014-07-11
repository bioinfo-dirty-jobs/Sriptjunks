setwd("/home/maurizio/Desktop/CountExperiments_Rnaseq/")

XA <-read.table("SampleA.count.txt")
Va <- read.table("hg19_toMaurizio_kb.txt")
valori <-as.numeric(as.character(XA[,1]))
totread=sum(valori)
totuniq=length(valori)

tofor (i in valori)
  re<-c(i/totuniq)