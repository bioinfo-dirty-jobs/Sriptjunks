xx<-as.character(az[,1])
yy<-as.character(az[,2])
match(xx,yy)
a<-read.table("dati_down.csv",sep="\t",header=T)