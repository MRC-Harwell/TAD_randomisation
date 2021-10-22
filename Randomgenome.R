#!/usr/bin/env Rscript

#Randomises the gene identity within each chromosome

for(i in seq(1,100)){
  genes=read.table("mm10_proteincodinggenes_biomart_sorted.bed")

  
  genes$V5 <- ave(genes$V5, genes$V1, FUN = sample)
  
  writefilestring=paste("mm10_RANDOMISEDgenome",i,".bed",sep="")
  write.table(genes,file=writefilestring,quote=F,sep="\t" ,row.names=F,col.names=F)
  
}

