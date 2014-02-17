Args<-commandArgs()
featureX<-Args[5]

x=read.table(featureX)
#y=as.table(matrix(1,dim(x)[1],4))

library(CENTIPEDE)
centFit <- fitCentipede(Xlist = list(DNase=as.matrix(x+0.1)))
write.table(centFit$PostPr,file= paste(featureX,"prob",sep="."),quote = FALSE,row.names=FALSE,col.names=FALSE)
