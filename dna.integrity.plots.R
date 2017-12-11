#Code needs to do the following:
#notes: only look at smokers.
#using 10 Mut/Mb as cutoff
#

#Figure 1 - Landscape of alteration of genome stability related proteins in non-small cell lung cancer.


#code to make figure 1

load('/Users/michaelsharpnack/Desktop/Kai/nsclc.11202017.RData')
data.all <- cbind(rbind(patient.mutburden.luad,neoantigens.burden.luad,as.numeric(smoking.luad[,3]),mm.sig.path.luad+mm.sig.mis.path.luad),
                  rbind(patient.mutburden.lusc,neoantigens.burden.lusc,as.numeric(smoking.lusc[,3]),mm.sig.path.lusc+mm.sig.mis.path.lusc))
data.all[4:20,][data.all[4:20,] > 2] <- 2
data.all[3,data.all[3,] == 1] <- 0
data.all[3,data.all[3,] > 1] <- 1
data.all[1,] <- log2(data.all[1,])/log2(max(data.all[1,]))
data.all[2,] <- data.all[2,]/max(data.all[2,])
data.all.temp <- data.all
data.all.temp[4:20,][data.all.temp[4:20,] > 1] <- 1
data.all[,1:461] <- data.all[,order(colSums(data.all.temp[4:20,1:461]),decreasing=TRUE)]
data.all[,462:926] <- data.all[,order(colSums(data.all.temp[4:20,462:926]),decreasing=TRUE)+461]
data.all <- rbind(c(rep(0,461),rep(1,465)),data.all)
rownames(data.all) <- c('Histology','Mutation Burden','Neoantigen Burden','Smoking Status',unique(genes.sig[,1]))
data.all <- data.all[c(rownames(data.all)[1:4],names(sort(rowSums(data.all[5:21,]),decreasing=TRUE))),]
data.all[2,] <- data.all[2,]+2
data.all[3,] <- data.all[3,]+4
data.all[4,] <- data.all[4,]+6
data.all[5:21,] <- data.all[5:21,]+8

write.table(data.all[c(2,4:21),1:461],file='/Users/michaelsharpnack/Desktop/Kai/data.luad.tab',sep='\t')
write.table(data.all[c(2,4:21),462:926],file='/Users/michaelsharpnack/Desktop/Kai/data.lusc.tab',sep='\t')

#Supplementary Figure 1



#Figure 2 - Association between smoking and tumor mutation burden.
figure2.plot <- function(patient.mutburden,smoking){
  #plot histogram with smoker and nonsmokers together
  hist(patient.mutburden[smoking[,3] > 1],col='skyblue',border=F,50)
  hist(patient.mutburden[smoking[,3] == 1],add=T,col=scales::alpha('red',.5),border=F,10)
  t.test(patient.mutburden[smoking[,3] == 1],patient.mutburden[smoking[,3] > 1])

  color <- 'blue'
  #plot association between TMB & smoking
  plot(as.numeric(smoking[,1]),patient.mutburden,pch=19,col=color)
  abline(lm(patient.mutburden~as.numeric(smoking[,1])),lwd=2)
  cor.test(as.numeric(smoking[is.na(smoking[,1]) == FALSE,1]),patient.mutburden[is.na(smoking[,1]) == FALSE],method='spearman')
}

load('/Users/michaelsharpnack/Desktop/Kai/luad.10022017.RData')
figure2.plot(patient.mutburden,smoking)
mutburden.smok.luad <- patient.mutburden[smoking[,3] > 1]
mutburden.nonsmok.luad <- patient.mutburden[smoking[,3] == 1]
mutburden.luad <- patient.mutburden
load('/Users/michaelsharpnack/Desktop/Kai/lusc.10022017.RData')
figure2.plot(patient.mutburden,smoking)
mutburden.smok.lusc <- patient.mutburden[smoking[,3] > 1]
mutburden.nonsmok.lusc <- patient.mutburden[smoking[,3] == 1]
mutburden.lusc <- patient.mutburden
#There is no difference between mutation burden in LUSC and LUAD smokers
t.test(mutburden.smok.luad,mutburden.smok.lusc)
t.test(mutburden.nonsmok.luad,mutburden.nonsmok.lusc)
t.test(mutburden.luad,mutburden.lusc)
t.test(mutburden.nonsmok.luad,mutburden.smok.luad)
t.test(mutburden.nonsmok.lusc,mutburden.smok.lusc)

#Figure 3 - Enrichment of genomic stability related genes and pathway inactivation in hypermutant tumors

fisher.function <- function(patient.mutburden,smoking,mm.sig,path.names,genes.sig,cna){
  #cutoff is 10 Mut/Mb
  names.temp1 <- rownames(smoking)[patient.mutburden < 10 & is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE]
  names.temp2 <- rownames(smoking)[patient.mutburden >= 10 & is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE]
  #only use smokers
  names.temp1 <- intersect(names.temp1,rownames(smoking)[smoking[,3] > 1])
  names.temp2 <- intersect(names.temp2,rownames(smoking)[smoking[,3] > 1])
  #test for pathway enrichment
  fisher.smokingtmb <- matrix(0,nrow=length(path.names),ncol=3)
  rownames(fisher.smokingtmb) <- path.names
  for(i in 1:length(path.names)){
    if(sum(genes.sig[,1] == path.names[i]) > 1){
      temp <- fisher.test(rbind(c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) > 0),
                                sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) == 0)),
                              c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) > 0),
                                sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) == 0))))
    } else {
      temp <- fisher.test(rbind(c(sum(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1 | cna[genes.sig[,1] == path.names[i],names.temp1] < -1),
                                  sum(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 0 & cna[genes.sig[,1] == path.names[i],names.temp1] >= -1)),
                                c(sum(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1 | cna[genes.sig[,1] == path.names[i],names.temp2] < -1),
                                  sum(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 0 & cna[genes.sig[,1] == path.names[i],names.temp2] >= -1))))
    }
    fisher.smokingtmb[i,1] <- temp$p.value
    fisher.smokingtmb[i,3] <- temp$estimate[1]
  }
  fisher.smokingtmb[,2] <- p.adjust(fisher.smokingtmb[,1],method='BH')
  #test for signature enrichment
  fisher.smoking.sigtmb <- matrix(0,nrow=dim(mm.sig)[1],ncol=3)
  rownames(fisher.smoking.sigtmb) <- rownames(mm.sig)
  for(i in 1:dim(mm.sig)[1]){
    temp <- fisher.test(rbind(c(sum(mm.sig[i,names.temp1] == 1 | cna[i,names.temp1] < -1),
                                sum(mm.sig[i,names.temp1] == 0 & cna[i,names.temp1] >= -1)),
                              c(sum(mm.sig[i,names.temp2] == 1 | cna[i,names.temp2] < -1),
                                sum(mm.sig[i,names.temp2] == 0 & cna[i,names.temp2] >= -1))))
    fisher.smoking.sigtmb[i,1] <- temp$p.value
    fisher.smoking.sigtmb[i,3] <- temp$estimate[1]
  }
  fisher.smoking.sigtmb[,2] <- p.adjust(fisher.smoking.sigtmb[,1],method='BH')
  
  plot(fisher.smoking.sigtmb[,3],-log10(fisher.smoking.sigtmb[,2]),xlim=c(0,2),ylim=c(0,2.8),pch=19,cex=1,col='blue')
  abline(v=1,h=1,lwd=2)
  
  plot(fisher.smokingtmb[,3],-log10(fisher.smokingtmb[,2]),xlim=c(0,2),ylim=c(0,2.8),pch=19,cex=2,col='blue')
  abline(v=1,h=1,lwd=2)

return(list("gene" = fisher.smoking.sigtmb,"pathway" = fisher.smokingtmb))
}

load('/Users/michaelsharpnack/Desktop/Kai/luad.10022017.RData')
fisher.matrix.luad <- fisher.function(patient.mutburden,smoking,mm.sig,path.names,genes.sig,cna)
data.luad <- list(patient.mutburden,smoking,mm.sig,path.names,genes.sig,cna)
load('/Users/michaelsharpnack/Desktop/Kai/lusc.10022017.RData')
fisher.matrix.lusc <- fisher.function(patient.mutburden,smoking,mm.sig,path.names,genes.sig,cna)
fisher.matrix.nsclc <- fisher.function(c(data.luad[[1]],patient.mutburden),rbind(data.luad[[2]],smoking),
                                      cbind(data.luad[[3]],mm.sig),path.names,genes.sig,cbind(data.luad[[6]],cna))

plot(-log10(fisher.matrix.luad$gene[,2]),-log10(fisher.matrix.lusc$gene[,2]),xlim=c(0,2.7),ylim=c(0,2.7),pch=19,col='blue')
lines(c(-1,3),c(-1,3),lwd=2)
plot(-log10(fisher.matrix.luad$pathway[,2]),-log10(fisher.matrix.lusc$pathway[,2]),xlim=c(0,2.7),ylim=c(0,2.7),pch=19,col='blue')
lines(c(-1,3),c(-1,3),lwd=2)

#plot the tumor burden of MSS and MSI-H tumors

load('/Users/michaelsharpnack/Desktop/Kai/nsclc.11202017.RData')
msi.stat <- fread('/Users/michaelsharpnack/Downloads/nm.4191-S3.csv')
msi.stat <- as.matrix(msi.stat)
rownames(msi.stat) <- paste(msi.stat[,1],'-01',sep='')
msi.stat.luad <- msi.stat[intersect(rownames(msi.stat),names(patient.mutburden.luad)),]
msi.stat.lusc <- msi.stat[intersect(rownames(msi.stat),names(patient.mutburden.lusc)),]
colnames(mm.sig.path.luad) <- names(patient.mutburden.luad)
colnames(mm.sig.path.lusc) <- names(patient.mutburden.lusc)
colnames(mm.sig.mis.path.luad) <- names(patient.mutburden.luad)
colnames(mm.sig.mis.path.lusc) <- names(patient.mutburden.lusc)

boxplot(patient.mutburden.luad[rownames(msi.stat.luad)[msi.stat.luad[,7] == 'MSS']],
        patient.mutburden.luad[rownames(msi.stat.luad)[msi.stat.luad[,7] != 'MSS']],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.path.luad[10,rownames(msi.stat.luad)] == 0]],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.path.luad[10,rownames(msi.stat.luad)] > 0]],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.mis.path.luad[10,rownames(msi.stat.luad)] == 0]],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.mis.path.luad[10,rownames(msi.stat.luad)] > 0]],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.path.luad[10,rownames(msi.stat.luad)] > 0 & msi.stat.luad[,7] == 'MSS']],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.path.luad[10,rownames(msi.stat.luad)] > 0 & msi.stat.luad[,7] != 'MSS']],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.path.luad[10,rownames(msi.stat.luad)] == 0 & msi.stat.luad[,7] == 'MSS']],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.path.luad[10,rownames(msi.stat.luad)] == 0 & msi.stat.luad[,7] != 'MSS']],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.mis.path.luad[10,rownames(msi.stat.luad)] == 0 & msi.stat.luad[,7] != 'MSS']],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.mis.path.luad[10,rownames(msi.stat.luad)] == 0 & msi.stat.luad[,7] != 'MSS']],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.mis.path.luad[10,rownames(msi.stat.luad)] == 0 & msi.stat.luad[,7] != 'MSS']],
        patient.mutburden.luad[rownames(msi.stat.luad)[mm.sig.mis.path.luad[10,rownames(msi.stat.luad)] == 0 & msi.stat.luad[,7] != 'MSS']]
        )

boxplot(patient.mutburden.lusc[rownames(msi.stat.lusc)[msi.stat.lusc[,7] == 'MSS']],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[msi.stat.lusc[,7] != 'MSS']],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.path.lusc[10,rownames(msi.stat.lusc)] == 0]],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.path.lusc[10,rownames(msi.stat.lusc)] > 0]],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.mis.path.lusc[10,rownames(msi.stat.lusc)] == 0]],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.mis.path.lusc[10,rownames(msi.stat.lusc)] > 0]],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.path.lusc[10,rownames(msi.stat.lusc)] > 0 & msi.stat.lusc[,7] == 'MSS']],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.path.lusc[10,rownames(msi.stat.lusc)] > 0 & msi.stat.lusc[,7] != 'MSS']],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.path.lusc[10,rownames(msi.stat.lusc)] == 0 & msi.stat.lusc[,7] == 'MSS']],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.path.lusc[10,rownames(msi.stat.lusc)] == 0 & msi.stat.lusc[,7] != 'MSS']],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.mis.path.lusc[10,rownames(msi.stat.lusc)] == 0 & msi.stat.lusc[,7] != 'MSS']],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.mis.path.lusc[10,rownames(msi.stat.lusc)] == 0 & msi.stat.lusc[,7] != 'MSS']],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.mis.path.lusc[10,rownames(msi.stat.lusc)] == 0 & msi.stat.lusc[,7] != 'MSS']],
        patient.mutburden.lusc[rownames(msi.stat.lusc)[mm.sig.mis.path.lusc[10,rownames(msi.stat.lusc)] == 0 & msi.stat.lusc[,7] != 'MSS']]
)











#Extras
#plot association between Neoantigen burden and TMB
plot(log2(patient.mutburden),neoantigens.burden,pch=19,col=color)
abline(lm(neoantigens.burden~log2(patient.mutburden)),lwd=2)
cor(log2(patient.mutburden),neoantigens.burden,method='spearman')

#plot association between Neoantigen burden and smoking pack years
plot(as.numeric(smoking.temp),neoantigens.burden,pch=19,col=color)
abline(lm(neoantigens.burden~as.numeric(smoking.temp)),lwd=2)
cor(as.numeric(smoking.temp[is.na(smoking.temp) == FALSE]),neoantigens.burden[is.na(smoking.temp) == FALSE],method='spearman')


