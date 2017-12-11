 ########################################################################################################################################################################  
#main function
library(data.table)
#2 == LUAD, 3 == LUSC
u = 2
#0 == no clonality info, 1 == clonality info
q = 0

tcga.maf.loader <- function(q,u){
  ##initialize variables##
  #load in the genes in the DNA integrity signature
  genes.sig <- fread('/Users/michaelsharpnack/Desktop/Kai/DNA%20repair%20related%20genes%20and%20DNA%20polymerases%20genes.csv')
  genes.symbol <- genes.sig[[2]]
  genes.sig <- as.matrix(genes.sig)
  rownames(genes.sig) <- genes.sig[,2]
  
  #do the MSigDb50 & purity/immune ttests
  immune.celltypes <- fread('/Users/michaelsharpnack/Downloads/immuneEstimation.txt')
  immune.celltypes <- as.matrix(immune.celltypes)
  rownames(immune.celltypes) <- immune.celltypes[,1]
  immune.celltypes <- immune.celltypes[,-1]
  class(immune.celltypes) <- 'numeric'
  
  ssgsea <- fread("/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/TCGA_Hallmark_gene_sets_ssGSEA_scores/tcga_fpkm_hallmarks_ssgsea.txt")
  ssgsea.rows <- ssgsea[[1]]
  ssgsea <- ssgsea[,-1]
  ssgsea <- as.matrix(ssgsea)
  class(ssgsea) <- 'numeric'
  ssgsea.samples <- read.csv("/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/TCGA_Hallmark_gene_sets_ssGSEA_scores/tcga_fpkm_samples.txt",stringsAsFactors=FALSE,sep='\t',header=FALSE)
  colnames(ssgsea) <- substr(ssgsea.samples[,1],1,15)
  rownames(ssgsea) <- ssgsea.rows
  rm(ssgsea.rows,ssgsea.samples)
  
  tumor.purity <- fread('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/ncomms9971-s2.csv')
  tumor.purity <- as.matrix(tumor.purity)
  rownames(tumor.purity) <- tumor.purity[,1]
  rownames(tumor.purity) <- substr(rownames(tumor.purity),1,15)
  tumor.purity <- tumor.purity[,-c(1,2,8)]
  class(tumor.purity) <- 'numeric'
  
  neoantigens <- fread('/Users/michaelsharpnack/Desktop/Kai/neoantigen_burden/nsclc.neoepitopes.csv')
  neoantigens <- as.matrix(neoantigens)
  rownames(neoantigens) <- neoantigens[,1]
  neoantigens <- neoantigens[,-1]
  neoantigens.burden <- as.numeric(neoantigens[,6])
  names(neoantigens.burden) <- paste(rownames(neoantigens),'-01',sep='')
  neoantigens.burden <- neoantigens.burden[is.na(neoantigens.burden) == FALSE]
  rm(neoantigens)
  
  #load in clinical data for tcga
  if(u == 2){
    clinical <- fread('/Users/michaelsharpnack/Downloads/gdac.broadinstitute.org_LUAD.Merge_Clinical.Level_1.2016012800.0.0/LUAD.clin.merged.txt')
  } else if(u == 3){
    clinical <- fread('/Users/michaelsharpnack/Downloads/gdac.broadinstitute.org_LUSC.Merge_Clinical.Level_1.2016012800.0.0/LUSC.clin.merged.txt')
  }
  clinical <- as.matrix(clinical)
  clinical <- t(clinical)
  colnames(clinical) <- clinical[1,]
  clinical <- clinical[-1,]
  smoking <- clinical[,grep('smok',colnames(clinical))]
  rownames(smoking) <- paste(clinical[,'patient.bcr_patient_barcode'],'-01',sep='')
  smoking <- smoking[is.na(smoking[,3]) == FALSE,]
  rownames(smoking) <- toupper(rownames(smoking))
  
  #load in CNV data
  if(u == 2){
    cancer <- 'luad'
  } else if(u == 3){
    cancer <- 'lusc'
  }
  load('/Users/michaelsharpnack/Desktop/Il-Jin/TP53/gtf.data.RData')
  temp.dir <- paste("/Users/michaelsharpnack/Desktop/Il-Jin/TP53/",cancer,"/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/",sep="")
  rnaseq.file <- dir(temp.dir)
  file.sample.map <- read.csv(paste("/Users/michaelsharpnack/Desktop/Il-Jin/TP53/",cancer,"/FILE_SAMPLE_MAP.txt",sep=""),stringsAsFactors=FALSE,header=FALSE,skip=1,sep="\t")
  file.sample.map.match <- substr(file.sample.map$V2,1,28)
  names(file.sample.map.match) <- file.sample.map$V1
  rna.file.barcode <- file.sample.map.match[intersect(names(file.sample.map.match),rnaseq.file)]
  patients <- substr(rna.file.barcode,1,15)
  patients = patients[as.numeric(substr(patients,14,15)) < 10]
  cnv <- matrix(0,nrow=20247,ncol=length(patients))
  rownames(cnv) <- rownames(gtf.temp.5)
  colnames(cnv) <- patients
  for(i in 1:length(patients)){
    cnv.temp <- read.csv(paste("/Users/michaelsharpnack/Desktop/Il-Jin/TP53/",cancer,"/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/",names(patients)[i],sep=""),sep="\t",stringsAsFactors=FALSE)
    cnv.temp <- cnv.temp[order(cnv.temp$Chromosome),]
    for(j in 1:dim(cnv.temp)[1]){
      if(nchar(cnv.temp$Start[j]) == 1){
        cnv.temp$Start[j] <- paste("00000000",cnv.temp$Start[j],sep="")
      }
      if(nchar(cnv.temp$Start[j]) == 2){
        cnv.temp$Start[j] <- paste("0000000",cnv.temp$Start[j],sep="")
      }
      if(nchar(cnv.temp$Start[j]) == 3){
        cnv.temp$Start[j] <- paste("000000",cnv.temp$Start[j],sep="")
      }
      if(nchar(cnv.temp$Start[j]) == 4){
        cnv.temp$Start[j] <- paste("00000",cnv.temp$Start[j],sep="")
      }
      if(nchar(cnv.temp$Start[j]) == 5){
        cnv.temp$Start[j] <- paste("0000",cnv.temp$Start[j],sep="")
      }
      if(nchar(cnv.temp$Start[j]) == 6){
        cnv.temp$Start[j] <- paste("000",cnv.temp$Start[j],sep="")
      }
      if(nchar(cnv.temp$Start[j]) == 7){
        cnv.temp$Start[j] <- paste("00",cnv.temp$Start[j],sep="")
      }
      if(nchar(cnv.temp$Start[j]) == 8){
        cnv.temp$Start[j] <- paste("0",cnv.temp$Start[j],sep="")
      }
      if(nchar(cnv.temp$End[j]) == 1){
        cnv.temp$End[j] <- paste("00000000",cnv.temp$End[j],sep="")
      }
      if(nchar(cnv.temp$End[j]) == 2){
        cnv.temp$End[j] <- paste("0000000",cnv.temp$End[j],sep="")
      }
      if(nchar(cnv.temp$End[j]) == 3){
        cnv.temp$End[j] <- paste("000000",cnv.temp$End[j],sep="")
      }
      if(nchar(cnv.temp$End[j]) == 4){
        cnv.temp$End[j] <- paste("00000",cnv.temp$End[j],sep="")
      }
      if(nchar(cnv.temp$End[j]) == 5){
        cnv.temp$End[j] <- paste("0000",cnv.temp$End[j],sep="")
      }
      if(nchar(cnv.temp$End[j]) == 6){
        cnv.temp$End[j] <- paste("000",cnv.temp$End[j],sep="")
      }
      if(nchar(cnv.temp$End[j]) == 7){
        cnv.temp$End[j] <- paste("00",cnv.temp$End[j],sep="")
      }
      if(nchar(cnv.temp$End[j]) == 8){
        cnv.temp$End[j] <- paste("0",cnv.temp$End[j],sep="")
      }
    }
    
    
    for(j in 1:dim(cnv.temp)[1]){
      cnv[rownames(gtf.temp.4[[which(unique(cnv.temp$Chromosome) == cnv.temp$Chromosome[j])]])[
        (gtf.temp.4[[which(unique(cnv.temp$Chromosome) == cnv.temp$Chromosome[j])]][,2] > cnv.temp$Start[j]) & 
          (gtf.temp.4[[which(unique(cnv.temp$Chromosome) == cnv.temp$Chromosome[j])]][,3] < cnv.temp$End[j])],i] <- 
        cnv.temp$Segment_Mean[j]
      
    }
    if(i %% 20 == 0){
      print(i)
    }
  }
  cnv <- cnv[,patients]
  
  break.points <- matrix(0,nrow = length(patients),ncol = 3)
  rownames(break.points) <- patients
  percent.altered <- matrix(0,nrow = length(patients),ncol = 3)
  rownames(percent.altered) <- patients
  for(i in 1:length(patients)){
    cnv.temp <- read.csv(paste("Desktop/Il-Jin/TP53/",cancer,"/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/",names(patients)[i],sep=""),sep="\t",stringsAsFactors=FALSE)
    
    break.points[i,1] <- sum(abs(cnv.temp$Segment_Mean) > 0.3)
    break.points[i,2] <- sum(cnv.temp$Segment_Mean < -0.3)
    break.points[i,3] <- sum(cnv.temp$Segment_Mean > 0.3)
    
    percent.altered[i,1] <- sum((cnv.temp$End-cnv.temp$Start)[abs(cnv.temp$Segment_Mean) > 0.3])/sum(as.numeric(cnv.temp$End-cnv.temp$Start))
    percent.altered[i,2] <- sum((cnv.temp$End-cnv.temp$Start)[cnv.temp$Segment_Mean < -0.3])/sum(as.numeric(cnv.temp$End-cnv.temp$Start))
    percent.altered[i,3] <- sum((cnv.temp$End-cnv.temp$Start)[cnv.temp$Segment_Mean > 0.3])/sum(as.numeric(cnv.temp$End-cnv.temp$Start))
  }
  break.points <- break.points[patients,]
  percent.altered <- percent.altered[patients,]
  rm(cnv.temp,file.sample.map,gtf,gtf.temp.2,gtf.temp.3,gtf.temp.5,cancer,file.sample.map.match,gtf.temp.4,i,j,patients,rna.file.barcode,rnaseq.file,temp.dir,gtf.read)
  
  #load in RNAseq data
  cancers <- c('BLCA','BRCA','CESC','COAD','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC')
  setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/')
  temp.rna.geneid <- read.csv(paste('tcga.rnaseq/',cancers[1],'.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',sep=""),stringsAsFactors=FALSE,sep='\t')[,1]
  temp.names <- strsplit(temp.rna.geneid,"\\|")
  temp.names <- temp.names[-c(1:30)]
  temp.names.2 <- vector(mode='numeric',length(temp.names))
  for(j in 1:length(temp.names)){
    temp.names.2[j] <- temp.names[[j]][1]
  }
  rna.geneid <- temp.names.2
  rm(temp.names,temp.names.2,temp.rna.geneid)
  setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query')
  if(u == 2){
    i = 12
  } else if(u == 3){
    i = 13
  }
  print(paste("Reading in RNAseq file for ",cancers[i]))
  rna <- fread(paste('tcga.rnaseq/',cancers[i],'.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',sep=""))
  rna <- as.matrix(rna[,-1])
  rna <- rna[-c(1:30),]
  class(rna) <- 'numeric'
  rownames(rna) <- rna.geneid
  colnames(rna) <- substr(colnames(rna),1,15)
  rm(cancers,rna.geneid)
  

  #load in the tcga maf files if they are individual
  #setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/tcga.maf/ovca')
  #maf.files <- dir()
  #for(i in 1:length(maf.files)){
  #  if(i == 1){
  #    temp <- fread(maf.files[i])
  #    mut.maf <- data.frame('Tumor_Sample_Barcode' = temp$Tumor_Sample_Barcode,'Hugo_Symbol' = temp$Hugo_Symbol,'Variant_Classification' = temp$Variant_Classification)
  #  } else {
  #    temp <- fread(maf.files[i])
  #    temp <- data.frame('Tumor_Sample_Barcode' = temp$Tumor_Sample_Barcode,'Hugo_Symbol' = temp$Hugo_Symbol,'Variant_Classification' = temp$Variant_Classification)
  #    mut.maf <- rbind(mut.maf,temp)
  #  }
  #}
  #colnames(mut.maf) <- c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification')
  
  #load in the tcga maf file
  setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/tcga.maf/')
  maf.files <- dir()
  mut.maf <- fread(maf.files[u])
  rm(maf.files)
  
  #mutation maf files that contain the clonality data
  if(q == 1){
    if(u == 2){
      luad.maf <- fread('/Users/michaelsharpnack/Downloads/tcga_pancancer_dcc_mafs_082115/mafs/tcga_luad_from_dcc.maf',select=c(1,9,11,12,13,16,42,80,81),skip=195093)
      luad.maf2 <- fread('/Users/michaelsharpnack/Downloads/tcga_pancancer_dcc_mafs_082115/mafs/tcga_luad_from_dcc.maf',select=c(1,9,11,12,13,16,42,80,81),nrows=195091)
      luad.maf.3 <- rbind(as.matrix(luad.maf2),as.matrix(luad.maf))
      mut.maf <- luad.maf.3
      rm(luad.maf2,luad.maf.3,luad.maf)
      mut.maf[,6] <- substr(mut.maf[,6],1,15)
      mut.maf <- cbind(mut.maf,as.numeric(mut.maf[,8])/(as.numeric(mut.maf[,8])+as.numeric(mut.maf[,9])))
      mut.maf <- cbind(mut.maf,vector(mode='numeric',dim(mut.maf)[1]))
      mut.maf <- cbind(mut.maf,vector(mode='numeric',dim(mut.maf)[1]))
      for(i in 1:dim(mut.maf)[1]){
        if(length(which(rownames(tumor.purity) == mut.maf[i,6])) > 0 & length(which(rownames(cnv) == mut.maf[i,1])) > 0){
          mut.maf[i,10] <- as.numeric(mut.maf[i,8])/(as.numeric(mut.maf[i,9])-as.numeric(mut.maf[i,9])*(1-tumor.purity[mut.maf[i,6],5])+as.numeric(mut.maf[i,8]))
          mut.maf[i,11] <- as.numeric(mut.maf[i,8])/(as.numeric(mut.maf[i,8])+as.numeric(mut.maf[i,9])/(10^cnv[mut.maf[i,1],mut.maf[i,6]]+1)*((10^cnv[mut.maf[i,1],mut.maf[i,6]]+1-2*(1-tumor.purity[mut.maf[i,6],5]))/tumor.purity[mut.maf[i,6],5]-2*(1-tumor.purity[mut.maf[i,6],5])^2))
        }
        if(i %% 10000 == 0){print(i)}
      }
    } else if(u == 3){
      lusc.maf <- fread('/Users/michaelsharpnack/Downloads/tcga_pancancer_dcc_mafs_082115/mafs/tcga_lusc_from_dcc.maf',select=c(1,9,11,12,13,16,80,81))
    }
  }
  mut.maf <- data.frame(mut.maf)
  
  #mut.types <- c('N Mutated','Total',unique(mut.maf$Variant_Classification),
  #               'A -> C','A -> G','A -> T','C -> A','C -> G','C -> T','G -> A','G -> C','G -> T','T -> A','T -> C','T -> G') 
  genes.symbol.int <- intersect(genes.symbol,unique(mut.maf$Hugo_Symbol))
  genes.sig <- genes.sig[genes.symbol.int,]
  path.names <- unique(genes.sig[,1])
  
  #coordinates for non-synonymous mutations in maf file
  #functional.muts <- union(union(union(union(which(mut.maf$Variant_Classification == 'Missense_Mutation'),which(mut.maf$Variant_Classification == 'Frame_Shift_Del')),
  #                                     union(which(mut.maf$Variant_Classification == 'Nonsense_Mutation'),which(mut.maf$Variant_Classification == 'Splice_Site'))),
  #                               union(which(mut.maf$Variant_Classification == 'Frame_Shift_Ins'),which(mut.maf$Variant_Classification == 'In_Frame_Del'))),
  #                       union(which(mut.maf$Variant_Classification == 'In_Frame_Ins'),which(mut.maf$Variant_Classification == 'Nonstop_Mutation')))
  functional.muts <- union(union(union(union(which(mut.maf$Variant_Classification == 'Frame_Shift_Del'),which(mut.maf$Variant_Classification == 'Nonsense_Mutation')),
                                 which(mut.maf$Variant_Classification == 'Frame_Shift_Ins')),
                                 which(mut.maf$Variant_Classification == 'In_Frame_Del')),
                           which(mut.maf$Variant_Classification == 'In_Frame_Ins'))
  functional.muts.mis <- union(union(which(mut.maf$Variant_Classification == 'Missense_Mutation'),which(mut.maf$Variant_Classification == 'Splice_Site')),
                                     which(mut.maf$Variant_Classification == 'Nonstop_Mutation'))
  
  #create mutation matrix for all genes
  #mm <- matrix(0,nrow=length(unique(mut.maf$Hugo_Symbol)),ncol=length(unique(mut.maf$Tumor_Sample_Barcode)))
  #genes.all <- unique(mut.maf$Hugo_Symbol)
  #rownames(mm) <- genes.all
  #colnames(mm) <- unique(mut.maf$Tumor_Sample_Barcode)
  #for(i in 1:length(genes.all)){
  #  mm[i,unique(mut.maf$Tumor_Sample_Barcode[intersect(which(mut.maf$Hugo_Symbol == genes.all[i]),functional.muts)])] <- 1
  #  if(i %% 1000 == 0){print(i)}
  #}
  #colnames(mm) <- substr(colnames(mm),1,15)
  #class(mm) <- 'numeric'
  
  #create mutation matrix for signature genes
  mm.sig <- matrix(0,nrow=length(genes.symbol.int),ncol=length(unique(mut.maf$Tumor_Sample_Barcode)))
  rownames(mm.sig) <- genes.symbol.int
  colnames(mm.sig) <- unique(mut.maf$Tumor_Sample_Barcode)
  for(i in 1:length(genes.symbol.int)){
    mm.sig[i,unique(mut.maf$Tumor_Sample_Barcode[intersect(which(mut.maf$Hugo_Symbol == genes.symbol.int[i]),functional.muts)])] <- 1
  }
  colnames(mm.sig) <- substr(colnames(mm.sig),1,15)
  
  
  mm.sig.mis <- matrix(0,nrow=length(genes.symbol.int),ncol=length(unique(mut.maf$Tumor_Sample_Barcode)))
  rownames(mm.sig.mis) <- genes.symbol.int
  colnames(mm.sig.mis) <- unique(mut.maf$Tumor_Sample_Barcode)
  for(i in 1:length(genes.symbol.int)){
    mm.sig.mis[i,unique(mut.maf$Tumor_Sample_Barcode[intersect(which(mut.maf$Hugo_Symbol == genes.symbol.int[i]),functional.muts.mis)])] <- 1
  }
  colnames(mm.sig.mis) <- substr(colnames(mm.sig.mis),1,15)
  
  #load in copy number alterations from cbioportal
  cna.files <- dir('/Users/michaelsharpnack/Desktop/Kai/gistic_cna')
  cna <- fread(paste('/Users/michaelsharpnack/Desktop/Kai/gistic_cna/',cna.files[u],sep=''))
  #cna <- rbind(cna,fread('/Users/michaelsharpnack/Desktop/Kai/luad_tcga_luad_tcga_gistic.2.txt'))
  cna.genes <- cna$GENE_ID
  cna <- as.matrix(cna[,3:dim(cna)[2]])
  rownames(cna) <- cna.genes
  cna <- cna[,colSums(is.nan(cna)) == 0]
  mm.sig <- mm.sig[,intersect(colnames(mm.sig),colnames(cna))]
  mm.sig.mis <- mm.sig.mis[,intersect(colnames(mm.sig),colnames(cna))]
  cna <- cna[,intersect(colnames(mm.sig),colnames(cna))]
  mm.sig <- mm.sig[intersect(rownames(mm.sig),rownames(cna)),]
  mm.sig.mis <- mm.sig.mis[intersect(rownames(mm.sig),rownames(cna)),]
  cna <- cna[intersect(rownames(mm.sig),rownames(cna)),]
  genes.sig <- genes.sig[intersect(rownames(mm.sig),rownames(cna)),]
  rm(cna.files,cna.genes,genes.symbol,genes.symbol.int)
  
  #calculate each tumor's mutation burden
  mut.maf$Tumor_Sample_Barcode <- substr(mut.maf$Tumor_Sample_Barcode,1,15)
  patient.mutburden <- vector(mode='numeric',length(unique(mut.maf$Tumor_Sample_Barcode)))
  for(j in 1:length(unique(mut.maf$Tumor_Sample_Barcode))){
    patient.mutburden[j] <- length(which(mut.maf$Tumor_Sample_Barcode == unique(mut.maf$Tumor_Sample_Barcode)[j]))
  }
  names(patient.mutburden) <- unique(mut.maf$Tumor_Sample_Barcode)
  patient.mutburden <- patient.mutburden[colnames(cna)]
  
  #load in Mb coverage for TCGA files
  beds.files <- dir('/Users/michaelsharpnack/Downloads/beds/')
  setwd('/Users/michaelsharpnack/Downloads/beds/')
  sample.basescovered <- vector(mode='numeric',length(beds.files))
  for(i in 1:length(beds.files)){
    beds.temp <- fread(beds.files[i])
    sample.basescovered[i] <- sum(as.numeric(beds.temp$V3)-as.numeric(beds.temp$V2))
  }
  names(sample.basescovered) <- beds.files
  
  beds.files <- dir('/Users/michaelsharpnack/Downloads/beds 2/')
  setwd('/Users/michaelsharpnack/Downloads/beds 2/')
  sample.basescovered.luad <- vector(mode='numeric',length(beds.files))
  for(i in 1:length(beds.files)){
    beds.temp <- fread(beds.files[i])
    sample.basescovered.luad[i] <- sum(as.numeric(beds.temp$V3)-as.numeric(beds.temp$V2))
  }
  names(sample.basescovered.luad) <- beds.files
  
  if(u == 2){
    patient.mutburden <- patient.mutburden/mean(sample.basescovered.luad)*1000000
  } else if(u == 3){
    patient.mutburden <- patient.mutburden/mean(sample.basescovered)*1000000
  }
  
  #remove the insane outlier in the lusc dataset
  if(u == 3){
    cna <- cna[,-which(names(patient.mutburden)==names(sort(patient.mutburden,decreasing=TRUE))[1])]
    mm.sig <- mm.sig[,-which(names(patient.mutburden)==names(sort(patient.mutburden,decreasing=TRUE))[1])]
    mm.sig.mis <- mm.sig.mis[,-which(names(patient.mutburden)==names(sort(patient.mutburden,decreasing=TRUE))[1])]
    patient.mutburden <- patient.mutburden[-which(names(patient.mutburden)==names(sort(patient.mutburden,decreasing=TRUE))[1])]
  }
  
  msi.stat <- fread('/Users/michaelsharpnack/Downloads/nm.4191-S3.csv')
  msi.stat <- as.matrix(msi.stat)
  rownames(msi.stat) <- paste(msi.stat[,1],'-01',sep='')
  
  
  #get the master patient list
  patients.temp <- intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(names(patient.mutburden),rownames(immune.celltypes)),rownames(smoking)),colnames(ssgsea)),rownames(tumor.purity)),colnames(mm.sig)),names(neoantigens.burden)),colnames(cnv)),colnames(rna))
  patient.mutburden <- patient.mutburden[patients.temp]
  immune.celltypes <- immune.celltypes[patients.temp,]
  smoking <- smoking[patients.temp,]
  ssgsea <- ssgsea[,patients.temp]
  tumor.purity <- tumor.purity[patients.temp,]
  mm.sig <- mm.sig[,patients.temp]
  mm.sig.mis <- mm.sig.mis[,patients.temp]
  neoantigens.burden <- neoantigens.burden[patients.temp]
  cnv <- cnv[,patients.temp]
  rna <- rna[,patients.temp]
  cna <- cna[,patients.temp]
  break.points <- break.points[patients.temp,]
  percent.altered <- percent.altered[patients.temp,]
  msi.stat <- msi.stat[intersect(rownames(msi.stat),patients.temp),]
  
  if(q == 1){
    #calculate each tumors clonal mutation burden
    clonal.mutburden <- vector(mode='numeric',length(unique(mut.maf[,6])))
    for(j in 1:length(unique(mut.maf[,6]))){
      clonal.mutburden[j] <- length(which(mut.maf[,6] == unique(mut.maf[,6])[j] & as.numeric(mut.maf[,10]) > 0.4))
    }
    names(clonal.mutburden) <- unique(mut.maf[,6])
    clonal.mutburden <- clonal.mutburden[patients.temp]
    #calculate each tumors subclonal mutation burden
    subclonal.mutburden <- vector(mode='numeric',length(unique(mut.maf[,6])))
    for(j in 1:length(unique(mut.maf[,6]))){
      subclonal.mutburden[j] <- length(which(mut.maf[,6] == unique(mut.maf[,6])[j] & as.numeric(mut.maf[,10]) < 0.4))
    }
    names(subclonal.mutburden) <- unique(mut.maf[,6])
    subclonal.mutburden <- subclonal.mutburden[patients.temp]
  }
  rm(q,u,i,j)
}

mm.sig.path <- matrix(0,nrow=17,ncol=dim(mm.sig)[2])
for(i in 1:17){
  if(length(which(genes.sig[,1] == unique(genes.sig[,1])[i])) > 1){
    mm.sig.path[i,] <- colSums(mm.sig[which(genes.sig[,1] == unique(genes.sig[,1])[i]),])
  }
  else{
    mm.sig.path[i,] <- mm.sig[which(genes.sig[,1] == unique(genes.sig[,1])[i]),]
  }
}

mm.sig.mis.path <- matrix(0,nrow=17,ncol=dim(mm.sig)[2])
for(i in 1:17){
  if(length(which(genes.sig[,1] == unique(genes.sig[,1])[i])) > 1){
    mm.sig.mis.path[i,] <- colSums(mm.sig.mis[which(genes.sig[,1] == unique(genes.sig[,1])[i]),])
  }
  else{
    mm.sig.mis.path[i,] <- mm.sig.mis[which(genes.sig[,1] == unique(genes.sig[,1])[i]),]
  }
}
mm.sig.mis.path[mm.sig.mis.path >= 1] <- 1

#################################################################################################################################################################
#test for association with DNA integrity inactivation and TMB/Neoantigens and TMB(NEO)/SPY

  #this is the original code used to produce the results used in the WCLC presentation
  #compute total mutation burden tumors with mutations in each signature gene
  mut.burden <- matrix(0,nrow=dim(cna)[1],ncol=6)
  rownames(mut.burden) <- rownames(cna)
  colnames(mut.burden) <- c('Mut+Del','Mut','Del','Missense','None in gene','Total')
  freq.altered <- mut.burden
  mut.pvalue <- mut.burden[,c(1:4)]
  for(i in 1:dim(cna)[1]){
    samples.mut <- colnames(mm.sig)[mm.sig[i,] == 1]
    samples.del <- colnames(cna)[cna[i,] < -1]
    samples.mis <- colnames(mm.sig.mis)[mm.sig.mis[i,] == 1]
  
    mut.burden[i,1] <- mean(patient.mutburden[union(samples.mut,samples.del)],na.rm=TRUE)
    mut.burden[i,2] <- mean(patient.mutburden[samples.mut],na.rm=TRUE)
    mut.burden[i,3] <- mean(patient.mutburden[samples.del],na.rm=TRUE)
    mut.burden[i,4] <- mean(patient.mutburden[samples.mis],na.rm=TRUE)
    mut.burden[i,5] <- mean(patient.mutburden[setdiff(names(patient.mutburden),union(samples.mut,samples.del))],na.rm=TRUE)

    freq.altered[i,1] <- length(union(samples.mut,samples.del))
    freq.altered[i,2] <- length(samples.mut)
    freq.altered[i,3] <- length(samples.del)
    freq.altered[i,4] <- length(samples.mis)
    freq.altered[i,5] <- length(setdiff(names(patient.mutburden),union(samples.mut,samples.del)))
    
    mut.pvalue[i,1] <- try(t.test(patient.mutburden[union(samples.mut,samples.del)],mu=mean(patient.mutburden))$p.value)
    mut.pvalue[i,2] <- try(t.test(patient.mutburden[samples.mut],mu=mean(patient.mutburden))$p.value)
    mut.pvalue[i,3] <- try(t.test(patient.mutburden[samples.del],mu=mean(patient.mutburden))$p.value)
    mut.pvalue[i,4] <- try(t.test(patient.mutburden[samples.mis],mu=mean(patient.mutburden))$p.value)

  }
  mut.burden[,6] <- mean(patient.mutburden)
  freq.altered[,6] <- length(patient.mutburden)
  mut.pvalue[nchar(mut.pvalue)> 100] <- NA
  
  #2-sided t-test between GENE inactivated patients and activated patients mutation and neoantigen burdens
  mut.burden.sig <- matrix(0,nrow=dim(cna)[1],ncol=9)
  rownames(mut.burden.sig) <- rownames(cna)
  for(i in 1:dim(cna)[1]){
    samples.mut <- colnames(mm.sig)[mm.sig[i,] == 1 | cna[i,] < -1]
    mut.burden.sig[i,1] <- length(samples.mut)
    mut.burden.sig[i,2] <- mean(patient.mutburden[samples.mut],na.rm=TRUE)
    mut.burden.sig[i,3] <- mean(patient.mutburden[setdiff(names(patient.mutburden),samples.mut)],na.rm=TRUE)
    mut.burden.sig[i,4] <- mean(neoantigens.burden[samples.mut],na.rm = TRUE)
    mut.burden.sig[i,5] <- mean(neoantigens.burden[setdiff(names(patient.mutburden),samples.mut)],na.rm = TRUE)
    mut.burden.sig[i,6] <- try(t.test(patient.mutburden[samples.mut],patient.mutburden[setdiff(names(patient.mutburden),samples.mut)])$p.value,silent = TRUE)
    mut.burden.sig[i,8] <- try(t.test(neoantigens.burden[samples.mut],neoantigens.burden[setdiff(names(neoantigens.burden),samples.mut)])$p.value,silent = TRUE)
  }
  class(mut.burden.sig) <- 'numeric'
  mut.burden.sig[,7] <- p.adjust(mut.burden.sig[,6],method='BH')
  mut.burden.sig[,9] <- p.adjust(mut.burden.sig[,8],method='BH')
  View(mut.burden.sig[is.na(mut.burden.sig[,6]) == FALSE,][order(mut.burden.sig[is.na(mut.burden.sig[,6]) == FALSE,8]),])
  
  #2-sided t-test between GENE inactivated patients and activated patients mutation and neoantigen burdens NORMALIZED BY SPY
  smoking.temp <- smoking[,1]
  smoking.temp[is.na(smoking[,1]) & smoking[,3] == 1] <- 0
  smoking.temp <- as.numeric(smoking.temp)
  names(smoking.temp) <- rownames(smoking)
  smoking.temp <- smoking.temp+1
  mut.burden.sig <- matrix(0,nrow=dim(cna)[1],ncol=9)
  rownames(mut.burden.sig) <- rownames(cna)
  for(i in 1:dim(cna)[1]){
    samples.mut <- colnames(mm.sig)[mm.sig[i,] == 1 | cna[i,] < -1]
    mut.burden.sig[i,1] <- length(samples.mut)
    mut.burden.sig[i,2] <- mean((patient.mutburden/smoking.temp)[samples.mut],na.rm=TRUE)
    mut.burden.sig[i,3] <- mean((patient.mutburden/smoking.temp)[setdiff(names(patient.mutburden),samples.mut)],na.rm=TRUE)
    mut.burden.sig[i,4] <- mean((neoantigens.burden/smoking.temp)[samples.mut],na.rm = TRUE)
    mut.burden.sig[i,5] <- mean((neoantigens.burden/smoking.temp)[setdiff(names(patient.mutburden),samples.mut)],na.rm = TRUE)
    mut.burden.sig[i,6] <- try(t.test((patient.mutburden/smoking.temp)[samples.mut][is.na((patient.mutburden/smoking.temp)[samples.mut]) == FALSE],(patient.mutburden/smoking.temp)[setdiff(names(patient.mutburden),samples.mut)][is.na((patient.mutburden/smoking.temp)[setdiff(names(patient.mutburden),samples.mut)]) == FALSE])$p.value,silent = TRUE)
    mut.burden.sig[i,8] <- try(t.test((neoantigens.burden/smoking.temp)[samples.mut][is.na((neoantigens.burden/smoking.temp)[samples.mut]) == FALSE],(neoantigens.burden/smoking.temp)[setdiff(names(neoantigens.burden),samples.mut)][is.na((neoantigens.burden/smoking.temp)[setdiff(names(neoantigens.burden),samples.mut)]) == FALSE])$p.value,silent = TRUE)
  }
  class(mut.burden.sig) <- 'numeric'
  mut.burden.sig[,7] <- p.adjust(mut.burden.sig[,6],method='BH')
  mut.burden.sig[,9] <- p.adjust(mut.burden.sig[,8],method='BH')
  View(mut.burden.sig[is.na(mut.burden.sig[,6]) == FALSE,][order(mut.burden.sig[is.na(mut.burden.sig[,6]) == FALSE,8]),])
  
  #calculate the mutation burden p-values for sig pathways
  mut.burden.path <- matrix(0,nrow=length(path.names),ncol=9)
  rownames(mut.burden.path) <- path.names
  for(i in 1:length(path.names)){
    if(length(rownames(genes.sig)[genes.sig[,1] == path.names[i]]) > 1){
      samples.mut <- colnames(mm.sig)[colSums(mm.sig[rownames(genes.sig)[genes.sig[,1] == path.names[i]],] == 1 | (cna[rownames(genes.sig)[genes.sig[,1] == path.names[i]],] < -1)) > 0]
    } else {
      samples.mut <- colnames(mm.sig)[mm.sig[rownames(genes.sig)[genes.sig[,1] == path.names[i]],] == 1 | cna[rownames(genes.sig)[genes.sig[,1] == path.names[i]],] < -1]
    }
    mut.burden.path[i,1] <- length(samples.mut)
    mut.burden.path[i,2] <- mean(patient.mutburden[samples.mut],na.rm=TRUE)
    mut.burden.path[i,3] <- mean(patient.mutburden[setdiff(names(patient.mutburden),samples.mut)],na.rm=TRUE)
    mut.burden.path[i,4] <- mean(neoantigens.burden[samples.mut],na.rm = TRUE)
    mut.burden.path[i,5] <- mean(neoantigens.burden[setdiff(names(patient.mutburden),samples.mut)],na.rm = TRUE)
    mut.burden.path[i,6] <- try(t.test(patient.mutburden[samples.mut],patient.mutburden[setdiff(names(patient.mutburden),samples.mut)])$p.value,silent = TRUE)
    mut.burden.path[i,8] <- try(t.test(neoantigens.burden[samples.mut],neoantigens.burden[setdiff(names(neoantigens.burden),samples.mut)])$p.value,silent = TRUE)
  }
  class(mut.burden.path) <- 'numeric'
  mut.burden.path[,7] <- p.adjust(mut.burden.path[,6],method='BH')
  mut.burden.path[,9] <- p.adjust(mut.burden.path[,8],method='BH')
  
  #calculate the mutation burden p-values for sig pathways / SPY
  mut.burden.path <- matrix(0,nrow=length(path.names),ncol=9)
  rownames(mut.burden.path) <- path.names
  for(i in 1:length(path.names)){
    if(length(rownames(genes.sig)[genes.sig[,1] == path.names[i]]) > 1){
      samples.mut <- colnames(mm.sig)[colSums(mm.sig[rownames(genes.sig)[genes.sig[,1] == path.names[i]],] == 1 | (cna[rownames(genes.sig)[genes.sig[,1] == path.names[i]],] < -1)) > 0]
    } else {
      samples.mut <- colnames(mm.sig)[mm.sig[rownames(genes.sig)[genes.sig[,1] == path.names[i]],] == 1 | cna[rownames(genes.sig)[genes.sig[,1] == path.names[i]],] < -1]
    }
    mut.burden.path[i,1] <- length(samples.mut)
    mut.burden.path[i,2] <- mean((patient.mutburden/smoking.temp)[samples.mut],na.rm=TRUE)
    mut.burden.path[i,3] <- mean((patient.mutburden/smoking.temp)[setdiff(names(patient.mutburden),samples.mut)],na.rm=TRUE)
    mut.burden.path[i,4] <- mean((neoantigens.burden/smoking.temp)[samples.mut],na.rm = TRUE)
    mut.burden.path[i,5] <- mean((neoantigens.burden/smoking.temp)[setdiff(names(patient.mutburden),samples.mut)],na.rm = TRUE)
    mut.burden.path[i,6] <- try(t.test((patient.mutburden/smoking.temp)[samples.mut][is.na((patient.mutburden/smoking.temp)[samples.mut]) == FALSE],(patient.mutburden/smoking.temp)[setdiff(names(patient.mutburden),samples.mut)][is.na((patient.mutburden/smoking.temp)[setdiff(names(patient.mutburden),samples.mut)]) == FALSE])$p.value,silent = TRUE)
    mut.burden.path[i,8] <- try(t.test((neoantigens.burden/smoking.temp)[samples.mut][is.na((neoantigens.burden/smoking.temp)[samples.mut]) == FALSE],(neoantigens.burden/smoking.temp)[setdiff(names(neoantigens.burden),samples.mut)][is.na((neoantigens.burden/smoking.temp)[setdiff(names(neoantigens.burden),samples.mut)]) == FALSE])$p.value,silent = TRUE)
  }
  class(mut.burden.path) <- 'numeric'
  mut.burden.path[,7] <- p.adjust(mut.burden.path[,6],method='BH')
  mut.burden.path[,9] <- p.adjust(mut.burden.path[,8],method='BH')
  
##################################################################################################################################
#model the mutation burden with DNA integrity genes as a penalized linear regression
library(glmnet)
#initialize the input matrix to the regression
input.matrix <- matrix(0,nrow=dim(cna)[2],ncol=dim(cna)[1])
input.matrix[t(cna) < -1] <- 1
input.matrix[t(mm.sig) > 0] <- 1
colnames(input.matrix) <- rownames(cna)
rownames(input.matrix) <- colnames(cna)

#k = number of cross validation for signature performance, holdout = number of samples in test set
k = 1000
holdout = 50
dnaintegrity.cv <- matrix(0,nrow=dim(input.matrix)[2],ncol=k)
rownames(dnaintegrity.cv) <- colnames(input.matrix)
cv.lasso.object <- list()
real.train <- matrix(0,nrow=dim(input.matrix)[1]-holdout,ncol=k)
predicted.train <- matrix(0,nrow=dim(input.matrix)[1]-holdout,ncol=k)
real.test <- matrix(0,nrow=holdout,ncol=k)
predicted.test <- matrix(0,nrow=holdout,ncol=k)
for(i in 1:k){
  sample.temp <- sample(1:dim(input.matrix)[1],dim(input.matrix)[1]-holdout)
  cv.lasso.object[[i]] <- cv.glmnet(input.matrix[sample.temp,],log2(patient.mutburden+1)[sample.temp])
  lasso.object <- glmnet(input.matrix[sample.temp,],log2(patient.mutburden+1)[sample.temp],family='poisson',lambda = cv.lasso.object[[i]]$lambda.min)
  dnaintegrity.cv[,i] <- lasso.object$beta[,1]
  real.train[,i] <- log2(patient.mutburden+1)[sample.temp]
  predicted.train[,i] <- predict(lasso.object,input.matrix[sample.temp,])
  real.test[,i] <- log2(patient.mutburden+1)[setdiff(1:dim(input.matrix)[1],sample.temp)]
  predicted.test[,i] <- predict(lasso.object,input.matrix[setdiff(1:dim(input.matrix)[1],sample.temp),])
  if(i %% 50 == 0){print(i)}
}

#checking to see how well the test/training predictions correlate with the actual mutation burdens
cor.train <- vector(mode='numeric',k)
cor.tester <- vector(mode='numeric',k)
for(i in 1:k){
  cor.tester[i] <- cor(predicted.test[,i],real.test[,i])
  cor.train[i] <- cor(real.train[,i],predicted.train[,i])
}

smoking.temp <- smoking[,1]
smoking.temp[is.na(smoking[,1]) & smoking[,3] == 1] <- 0
smoking.temp <- as.numeric(smoking.temp)
names(smoking.temp) <- rownames(smoking)
smoking.temp <- smoking.temp+1
k = 100
dnaintegrity.cv.spy <- matrix(0,nrow=150,ncol=k)
rownames(dnaintegrity.cv.spy) <- rownames(cna)
cv.lasso.object.spy <- list()
for(i in 1:k){

  cv.lasso.object.spy[[i]] <- cv.glmnet(input.matrix[is.na((patient.mutburden/smoking.temp)) == FALSE,],log2((patient.mutburden/smoking.temp)[is.na((patient.mutburden/smoking.temp)) == FALSE]+1))
  lasso.object <- glmnet(input.matrix[is.na((patient.mutburden/smoking.temp)) == FALSE,],log2((patient.mutburden/smoking.temp)[is.na((patient.mutburden/smoking.temp)) == FALSE]+1),family='poisson',lambda = cv.lasso.object.spy[[i]]$lambda.min)
  dnaintegrity.cv.spy[,i] <- lasso.object$beta[,1]
}



##################################################################################################################################
#make boxplots of significance for the association between gene inactivation and average TMB

luad.boxplot <- list()
luad.boxplot.names <- rownames(mm.sig)[mut.burden[,1] > 400][order(mut.pvalue[mut.burden[,1] > 400,1])][1:10]
for(i in 1:10){
  luad.boxplot[[i]] <- patient.mutburden[mm.sig[rownames(mm.sig)[mut.burden[,1] > 400][order(mut.pvalue[mut.burden[,1] > 400,1])][i],] == 1 | cna[rownames(mm.sig)[mut.burden[,1] > 400][order(mut.pvalue[mut.burden[,1] > 400,1])][i],] < -1]
}
boxplot(luad.boxplot,names = luad.boxplot.names,col=c(rep('skyblue',4),rep('white',6)))
abline(h=mean(patient.mutburden),lwd=2)

##################################################################################################################################
#plot histogram with smoker and nonsmokers together
hist(log2(patient.mutburden[smoking[,3] > 1]),col='skyblue',border=F,50)
hist(log2(patient.mutburden[smoking[,3] == 1]),add=T,col=scales::alpha('red',.5),border=F,30)
t.test(patient.mutburden[smoking[,3] == 1],patient.mutburden[smoking[,3] > 1])

color <- 'skyblue'
#plot association between TMB & smoking
smoking.temp <- smoking[,1]
names(smoking.temp) <- rownames(smoking)
smoking.temp[is.na(smoking[,1]) & smoking[,3] == 1] <- 0
plot(as.numeric(smoking.temp),log2(patient.mutburden),pch=19,col=color)
abline(lm(log2(patient.mutburden)~as.numeric(smoking.temp)),lwd=2)
cor.test(as.numeric(smoking.temp[is.na(smoking.temp) == FALSE]),log2(patient.mutburden[is.na(smoking.temp) == FALSE]),method='spearman')

#plot association between Neoantigen burden and TMB
plot(log2(patient.mutburden),neoantigens.burden,pch=19,col=color)
abline(lm(neoantigens.burden~log2(patient.mutburden)),lwd=2)
cor(log2(patient.mutburden),neoantigens.burden,method='spearman')

#plot association between Neoantigen burden and smoking pack years
plot(as.numeric(smoking.temp),neoantigens.burden,pch=19,col=color)
abline(lm(neoantigens.burden~as.numeric(smoking.temp)),lwd=2)
cor(as.numeric(smoking.temp[is.na(smoking.temp) == FALSE]),neoantigens.burden[is.na(smoking.temp) == FALSE],method='spearman')

#plot association between smoking pack years & tumor purity
plot(as.numeric(smoking.temp),tumor.purity[,4],pch=19,col='red')
cor(as.numeric(smoking.temp[is.na(smoking.temp) == FALSE & is.na(tumor.purity[,4]) == FALSE]),tumor.purity[,4][is.na(smoking.temp) == FALSE & is.na(tumor.purity[,4]) == FALSE],method='spearman')

#plot association between smoking pack years & immune cell subsets
cor.immunecelltypes.smokingpackyears <- vector(mode='numeric',6)
names(cor.immunecelltypes.smokingpackyears) <- colnames(immune.celltypes)
for(i in 1:6){
  cor.immunecelltypes.smokingpackyears[i] <- cor(as.numeric(smoking.temp[is.na(smoking.temp) == FALSE]),immune.celltypes[,i][is.na(smoking.temp) == FALSE],method='spearman')
}

#plot association between immune cell subsets and TMB
cor.immunecelltypes.tmb <- vector(mode='numeric',6)
names(cor.immunecelltypes.tmb) <- colnames(immune.celltypes)
for(i in 1:6){
  cor.immunecelltypes.tmb[i] <- cor(patient.mutburden,immune.celltypes[,i],method='spearman')
}
plot(log2(patient.mutburden),immune.celltypes[,3]/immune.celltypes[,2],method='spearman')


#correlation between RNA values and smoking pack years
#NOTE that none of these are significant after BH correction
cor.rna.smokingpackyears <- matrix(0,nrow=dim(rna)[1],ncol=3)
rownames(cor.rna.smokingpackyears) <- rownames(rna)
colnames(cor.rna.smokingpackyears) <- c('Rho','p-value','BH q-value')
for(i in 1:dim(rna)[1]){
  temp <- cor.test(as.numeric(smoking.temp[is.na(smoking.temp) == FALSE]),rna[i,][is.na(smoking.temp) == FALSE],method='spearman')
  cor.rna.smokingpackyears[i,1] <- temp$estimate
  cor.rna.smokingpackyears[i,2] <- temp$p.value
}
cor.rna.smokingpackyears[,3] <- p.adjust(cor.rna.smokingpackyears[,2],method='BH')
View(cor.rna.smokingpackyears[order(cor.rna.smokingpackyears[,2]),])

#correlation between RNA values and TMB
cor.rna.tmb <- matrix(0,nrow=dim(rna)[1],ncol=3)
rownames(cor.rna.tmb) <- rownames(rna)
colnames(cor.rna.tmb) <- c('Rho','p-value','BH q-value')
for(i in 1:dim(rna)[1]){
  temp <- cor.test(patient.mutburden,rna[i,],method='spearman')
  cor.rna.tmb[i,1] <- temp$estimate
  cor.rna.tmb[i,2] <- temp$p.value
}
cor.rna.tmb[,3] <- p.adjust(cor.rna.tmb[,2],method='BH')
View(cor.rna.tmb[order(cor.rna.tmb[,2]),])

#correlation between RNA values and NEO
cor.rna.neo <- matrix(0,nrow=dim(rna)[1],ncol=3)
rownames(cor.rna.neo) <- rownames(rna)
colnames(cor.rna.neo) <- c('Rho','p-value','BH q-value')
for(i in 1:dim(rna)[1]){
  temp <- cor.test(neoantigens.burden,rna[i,],method='spearman')
  cor.rna.neo[i,1] <- temp$estimate
  cor.rna.neo[i,2] <- temp$p.value
}
cor.rna.neo[,3] <- p.adjust(cor.rna.neo[,2],method='BH')
View(cor.rna.neo[order(cor.rna.neo[,2]),])
plot(cor.rna.neo[,1],-log10(cor.rna.neo[,3]))
hist(cor.rna.neo[,1],1000)

#correlation between RNA values and TMB in smokers vs nonsmokers
cor.rna.tmb.SmokerNonSmoker <- matrix(0,nrow=dim(rna)[1],ncol=6)
rownames(cor.rna.tmb.SmokerNonSmoker) <- rownames(rna)
colnames(cor.rna.tmb.SmokerNonSmoker) <- c('Smoker Rho','Smoker p-value','Smoker BH q-value','Non-Smoker Rho','Non-Smoker p-value','Non-Smoker BH q-value')
for(i in 1:dim(rna)[1]){
  temp <- cor.test(patient.mutburden[smoking[,3] > 1],rna[i,][smoking[,3] > 1],method='spearman')
  cor.rna.tmb.SmokerNonSmoker[i,1] <- temp$estimate
  cor.rna.tmb.SmokerNonSmoker[i,2] <- temp$p.value
  temp <- cor.test(patient.mutburden[smoking[,3] == 1],rna[i,][smoking[,3] == 1],method='spearman')
  cor.rna.tmb.SmokerNonSmoker[i,4] <- temp$estimate
  cor.rna.tmb.SmokerNonSmoker[i,5] <- temp$p.value
  
}
cor.rna.tmb.SmokerNonSmoker[,3] <- p.adjust(cor.rna.tmb.SmokerNonSmoker[,2],method='BH')
cor.rna.tmb.SmokerNonSmoker[,6] <- p.adjust(cor.rna.tmb.SmokerNonSmoker[,5],method='BH')
View(cor.rna.tmb.SmokerNonSmoker[order(cor.rna.tmb.SmokerNonSmoker[,2]),])
View(cor.rna.tmb.SmokerNonSmoker[order(cor.rna.tmb.SmokerNonSmoker[,5]),])
plot(cor.rna.tmb.SmokerNonSmoker[,1],cor.rna.tmb.SmokerNonSmoker[,4],pch=19,col='skyblue',cex=0.5)
abline(lm(cor.rna.tmb.SmokerNonSmoker[,4]~cor.rna.tmb.SmokerNonSmoker[,1]),lwd=2)
abline(v=0,h=0,lwd=2)

#correlation between RNA values and neoantigen burden in smokers vs nonsmokers
cor.rna.neo.SmokerNonSmoker <- matrix(0,nrow=dim(rna)[1],ncol=6)
rownames(cor.rna.neo.SmokerNonSmoker) <- rownames(rna)
colnames(cor.rna.neo.SmokerNonSmoker) <- c('Smoker Rho','Smoker p-value','Smoker BH q-value','Non-Smoker Rho','Non-Smoker p-value','Non-Smoker BH q-value')
for(i in 1:dim(rna)[1]){
  temp <- cor.test(neoantigens.burden[smoking[,3] > 1],rna[i,][smoking[,3] > 1],method='spearman')
  cor.rna.neo.SmokerNonSmoker[i,1] <- temp$estimate
  cor.rna.neo.SmokerNonSmoker[i,2] <- temp$p.value
  temp <- cor.test(neoantigens.burden[smoking[,3] == 1],rna[i,][smoking[,3] == 1],method='spearman')
  cor.rna.neo.SmokerNonSmoker[i,4] <- temp$estimate
  cor.rna.neo.SmokerNonSmoker[i,5] <- temp$p.value
  
}
cor.rna.neo.SmokerNonSmoker[,3] <- p.adjust(cor.rna.neo.SmokerNonSmoker[,2],method='BH')
cor.rna.neo.SmokerNonSmoker[,6] <- p.adjust(cor.rna.neo.SmokerNonSmoker[,5],method='BH')
View(cor.rna.neo.SmokerNonSmoker[order(cor.rna.neo.SmokerNonSmoker[,2]),])
View(cor.rna.neo.SmokerNonSmoker[order(cor.rna.neo.SmokerNonSmoker[,5]),])
plot(cor.rna.neo.SmokerNonSmoker[,1],cor.rna.neo.SmokerNonSmoker[,4],pch=19,col='skyblue',cex=0.5)
cor(cor.rna.neo.SmokerNonSmoker[,1][is.na(cor.rna.neo.SmokerNonSmoker[,1]) == FALSE & is.na(cor.rna.neo.SmokerNonSmoker[,4]) == FALSE],cor.rna.neo.SmokerNonSmoker[,4][is.na(cor.rna.neo.SmokerNonSmoker[,1]) == FALSE & is.na(cor.rna.neo.SmokerNonSmoker[,4]) == FALSE])
abline(lm(cor.rna.neo.SmokerNonSmoker[,4]~cor.rna.neo.SmokerNonSmoker[,1]),lwd=2)
abline(v=0,h=0,lwd=2)

#correlation between ssgsea values and TMB in smokers vs nonsmokers
cor.ssgsea.tmb.SmokerNonSmoker <- matrix(0,nrow=dim(ssgsea)[1],ncol=6)
rownames(cor.ssgsea.tmb.SmokerNonSmoker) <- rownames(ssgsea)
colnames(cor.ssgsea.tmb.SmokerNonSmoker) <- c('Smoker Rho','Smoker p-value','Smoker BH q-value','Non-Smoker Rho','Non-Smoker p-value','Non-Smoker BH q-value')
for(i in 1:dim(ssgsea)[1]){
  temp <- cor.test(patient.mutburden[smoking[,3] > 1],ssgsea[i,][smoking[,3] > 1],method='spearman')
  cor.ssgsea.tmb.SmokerNonSmoker[i,1] <- temp$estimate
  cor.ssgsea.tmb.SmokerNonSmoker[i,2] <- temp$p.value
  temp <- cor.test(patient.mutburden[smoking[,3] == 1],ssgsea[i,][smoking[,3] == 1],method='spearman')
  cor.ssgsea.tmb.SmokerNonSmoker[i,4] <- temp$estimate
  cor.ssgsea.tmb.SmokerNonSmoker[i,5] <- temp$p.value
  
}
cor.ssgsea.tmb.SmokerNonSmoker[,3] <- p.adjust(cor.ssgsea.tmb.SmokerNonSmoker[,2],method='BH')
cor.ssgsea.tmb.SmokerNonSmoker[,6] <- p.adjust(cor.ssgsea.tmb.SmokerNonSmoker[,5],method='BH')
View(cor.ssgsea.tmb.SmokerNonSmoker[order(cor.ssgsea.tmb.SmokerNonSmoker[,2]),])
View(cor.ssgsea.tmb.SmokerNonSmoker[order(cor.ssgsea.tmb.SmokerNonSmoker[,5]),])
plot(cor.ssgsea.tmb.SmokerNonSmoker[,1],cor.ssgsea.tmb.SmokerNonSmoker[,4],pch=19,col='skyblue',cex=1)

#correlation between ssgsea values and Neoantigen burden in smokers vs nonsmokers
cor.ssgsea.neo.SmokerNonSmoker <- matrix(0,nrow=dim(ssgsea)[1],ncol=6)
rownames(cor.ssgsea.neo.SmokerNonSmoker) <- rownames(ssgsea)
colnames(cor.ssgsea.neo.SmokerNonSmoker) <- c('Smoker Rho','Smoker p-value','Smoker BH q-value','Non-Smoker Rho','Non-Smoker p-value','Non-Smoker BH q-value')
for(i in 1:dim(ssgsea)[1]){
  temp <- cor.test(neoantigens.burden[smoking[,3] > 1],ssgsea[i,][smoking[,3] > 1],method='spearman')
  cor.ssgsea.neo.SmokerNonSmoker[i,1] <- temp$estimate
  cor.ssgsea.neo.SmokerNonSmoker[i,2] <- temp$p.value
  temp <- cor.test(neoantigens.burden[smoking[,3] == 1],ssgsea[i,][smoking[,3] == 1],method='spearman')
  cor.ssgsea.neo.SmokerNonSmoker[i,4] <- temp$estimate
  cor.ssgsea.neo.SmokerNonSmoker[i,5] <- temp$p.value
  
}
cor.ssgsea.neo.SmokerNonSmoker[,3] <- p.adjust(cor.ssgsea.neo.SmokerNonSmoker[,2],method='BH')
cor.ssgsea.neo.SmokerNonSmoker[,6] <- p.adjust(cor.ssgsea.neo.SmokerNonSmoker[,5],method='BH')
View(cor.ssgsea.neo.SmokerNonSmoker[order(cor.ssgsea.neo.SmokerNonSmoker[,2]),])
View(cor.ssgsea.neo.SmokerNonSmoker[order(cor.ssgsea.neo.SmokerNonSmoker[,5]),])
plot(cor.ssgsea.neo.SmokerNonSmoker[,1],cor.ssgsea.neo.SmokerNonSmoker[,4],pch=19,col='skyblue',cex=1)
abline(lm(cor.ssgsea.neo.SmokerNonSmoker[,4]~cor.ssgsea.neo.SmokerNonSmoker[,1]),lwd=2)
abline(v=0,h=0,lwd=2)



#correlation between RNA values and TMB/SPY
cor.rna.tmbspy <- matrix(0,nrow=dim(rna)[1],ncol=3)
rownames(cor.rna.tmbspy) <- rownames(rna)
colnames(cor.rna.tmbspy) <- c('Rho','p-value','BH q-value')
smoking.temp1 <- as.numeric(smoking.temp)+1
for(i in 1:dim(rna)[1]){
  temp <- cor.test((patient.mutburden/smoking.temp1)[is.na(smoking.temp1) == FALSE],rna[i,][is.na(smoking.temp1) == FALSE],method='spearman')
  cor.rna.tmbspy[i,1] <- temp$estimate
  cor.rna.tmbspy[i,2] <- temp$p.value
}
cor.rna.tmbspy[,3] <- p.adjust(cor.rna.tmbspy[,2],method='BH')
View(cor.rna.tmbspy[order(cor.rna.tmbspy[,2]),])


#copy number break points
#deleted vs. gained
plot(log2(break.points[,2]),log2(break.points[,3]),pch=19,col='red')
abline(lm(log2(break.points[,3])~log2(break.points[,2])),lwd=2)

#copy number and mutation burden
plot(log2(break.points[,1]),log2(patient.mutburden),pch=19,col='red')
abline(lm(log2(patient.mutburden)~log2(break.points[,1])),lwd=2)
cor(log2(patient.mutburden),log2(break.points[,1]),method='spearman')





##################################################################################################################################
#calculate ssgsea expression enrichment for each pathway in our signature
library(GSVA)
path.geneset <- list()
for(i in 1:length(path.names)){
  path.geneset[[i]] <- genes.sig[,2][genes.sig[,1] == path.names[i]][intersect(rownames(rna),genes.sig[,2][genes.sig[,1] == path.names[i]])]
}
ssgsea.integrity <- gsva(rna,path.geneset,method = 'ssgsea',rnaseq=TRUE)
rownames(ssgsea.integrity) <- path.names
ssgsea.integrity <- ssgsea.integrity[,patients.temp]


##################################################################################################################################

#make boxplot of immune cell concentrations in smokers/non-smokers and high/low TMB

names.temp1 <- rownames(smoking)[smoking[,3] == 1]
names.temp1 <- names.temp1[patient.mutburden[names.temp1] < 150]
names.temp2 <- rownames(smoking)[smoking[,3] == 1]
names.temp2 <- names.temp2[patient.mutburden[names.temp2] >= 150]

names.temp3 <- rownames(smoking)[smoking[,3] > 1]
names.temp3 <- names.temp3[patient.mutburden[names.temp3] < 150]
names.temp4 <- rownames(smoking)[smoking[,3] > 1]
names.temp4 <- names.temp4[patient.mutburden[names.temp4] >= 150]

names.temp5 <- rownames(smoking)[smoking[,3] == 1]
names.temp5 <- names.temp5[patient.mutburden[names.temp5] < 150]
names.temp6 <- rownames(smoking)[smoking[,3] > 1]
names.temp6 <- names.temp6[patient.mutburden[names.temp6] < 150]

names.temp7 <- rownames(smoking)[smoking[,3] == 1]
names.temp7 <- names.temp7[patient.mutburden[names.temp7] >= 150]
names.temp8 <- rownames(smoking)[smoking[,3] > 1]
names.temp8 <- names.temp8[patient.mutburden[names.temp8] >= 150]
for(i in 1:6){
  boxplot(immune.celltypes[names.temp1,i],immune.celltypes[names.temp2,i],immune.celltypes[names.temp3,i],immune.celltypes[names.temp4,i],immune.celltypes[names.temp5,i],immune.celltypes[names.temp6,i],immune.celltypes[names.temp7,i],immune.celltypes[names.temp8,i],
          at = c(1,2,4,5,7,8,10,11),main=colnames(immune.celltypes)[i])
  t.test(immune.celltypes[names.temp1,i],immune.celltypes[names.temp2,i])
}

##################################################################################################################################
#calculate gene lengths

hg19GeneLengths <- function(symbols)
{
  require(TxDb.Hsapiens.UCSC.hg19.knownGene) 
  require(org.Hs.eg.db)
  exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')    
  egs    = unlist(  mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG) )
  sapply(egs,function(eg)
  {
    exons = exons.db[[eg]]
    if(is.null(exons) == FALSE){
      exons = reduce(exons)
      sum( width(exons) )
    }
  })
}
gene.length <- hg19GeneLengths(rownames(mm))
gene.length.vec <- unlist(gene.length)

#exome length
gtf <- fread('/Users/michaelsharpnack/Desktop/Il-Jin/TP53/gtf/TCGA.hg19.June2011.gaf.txt')
gtf.exon <- gtf[which(gtf$FeatureType == 'gene'),]
gtf.exon <- gtf.exon[duplicated(gtf.exon$GeneLocus) == FALSE,]
gtf.exon <- gtf.exon[duplicated(gtf.exon$CompositeCoordinates) == FALSE,]
gene.length <- vector(mode='numeric',dim(gtf.exon)[1])
for(i in 1:dim(gtf.exon)[1]){
  gene.length[i] <- as.numeric(tail(strsplit(gtf.exon$FeatureCoordinates[i],"-")[[1]],n=1))
}
print(sum(gene.length))

##################################################################################################################################
#enrichment in TMB/SPY high vs. low for genes and pathways
names.temp1 <- rownames(smoking)[log2(patient.mutburden/as.numeric(smoking[,1])) < quantile(log2(patient.mutburden/as.numeric(smoking[,1]))[is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE])[3] & is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE]
names.temp2 <- rownames(smoking)[log2(patient.mutburden/as.numeric(smoking[,1])) >= quantile(log2(patient.mutburden/as.numeric(smoking[,1]))[is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE])[3] & is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE]
fisher.smoking <- matrix(0,nrow=length(path.names),ncol=3)
rownames(fisher.smoking) <- path.names
for(i in 1:length(path.names)){
  temp <- fisher.test(rbind(c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) > 0),
                              sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) == 0)),
                            c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) > 0),
                              sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) == 0))))
  print(path.names[i])
  print(temp)
  print(rbind(c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) > 0),
                sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) == 0)),
              c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) > 0),
                sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) == 0))))
  fisher.smoking[i,1] <- temp$p.value
  fisher.smoking[i,3] <- temp$estimate[1]
}
fisher.smoking[,2] <- p.adjust(fisher.smoking[,1],method='BH')

fisher.smoking.sig <- matrix(0,nrow=dim(mm.sig)[1],ncol=3)
rownames(fisher.smoking.sig) <- rownames(mm.sig)
for(i in 1:dim(mm.sig)[1]){
  temp <- fisher.test(rbind(c(sum(mm.sig[i,names.temp1] == 1 | cna[i,names.temp1] < -1),
                              sum(mm.sig[i,names.temp1] == 0 & cna[i,names.temp1] >= -1)),
                            c(sum(mm.sig[i,names.temp2] == 1 | cna[i,names.temp2] < -1),
                              sum(mm.sig[i,names.temp2] == 0 & cna[i,names.temp2] >= -1))))
  fisher.smoking.sig[i,1] <- temp$p.value
  fisher.smoking.sig[i,3] <- temp$estimate[1]
}
fisher.smoking.sig[,2] <- p.adjust(fisher.smoking.sig[,1],method='BH')

if(u == 2){
  plot(fisher.smoking[,3],-log10(fisher.smoking[,2]),xlim=c(0,2),ylim=c(0,5),pch=19,cex=2,col='skyblue')
  abline(v=1,h=1,lwd=2)
  plot(sort(log2(patient.mutburden/as.numeric(smoking[,1]))),col='skyblue',pch=19,cex=0.5,ylim=c(-3.5,10.5))
  abline(h=median(log2(patient.mutburden/as.numeric(smoking[,1])),na.rm=TRUE))
} else if(u == 3){
  plot(fisher.smoking[,3],-log10(fisher.smoking[,2]),xlim=c(0,2),ylim=c(0,5),pch=19,cex=2,col='red')
  abline(v=1,h=1,lwd=2)
  plot(sort(log2(patient.mutburden/as.numeric(smoking[,1]))),col='red',pch=19,cex=0.5,ylim=c(-3.5,10.5))
  abline(h=median(log2(patient.mutburden/as.numeric(smoking[,1])),na.rm=TRUE))
}
View(fisher.smoking[order(fisher.smoking[,1]),])

#fisher's test for enrichment of pathway inactivation in TMB high / TMB low

names.temp1 <- rownames(smoking)[patient.mutburden < quantile(patient.mutburden[is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE])[3] & is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE]
names.temp1 <- rownames(smoking)[patient.mutburden < 10 & is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE]

names.temp1 <- intersect(names.temp1,rownames(smoking)[smoking[,3] > 1])

names.temp2 <- rownames(smoking)[patient.mutburden >= quantile(patient.mutburden[is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE])[3] & is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE]
names.temp2 <- rownames(smoking)[patient.mutburden >= 10 & is.na(log2(patient.mutburden/as.numeric(smoking[,1]))) == FALSE]

names.temp2 <- intersect(names.temp2,rownames(smoking)[smoking[,3] > 1])

fisher.smokingtmb <- matrix(0,nrow=length(path.names),ncol=3)
rownames(fisher.smokingtmb) <- path.names
for(i in 1:length(path.names)){
  temp <- fisher.test(rbind(c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) > 0),
                              sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) == 0)),
                            c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) > 0),
                              sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) == 0))))
  print(path.names[i])
  print(temp)
  print(rbind(c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) > 0),
                sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp1] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp1] < -1) == 0)),
              c(sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) > 0 | colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) > 0),
                sum(colSums(mm.sig[genes.sig[,1] == path.names[i],names.temp2] == 1) == 0 & colSums(cna[genes.sig[,1] == path.names[i],names.temp2] < -1) == 0))))
  fisher.smokingtmb[i,1] <- temp$p.value
  fisher.smokingtmb[i,3] <- temp$estimate[1]
}
fisher.smokingtmb[,2] <- p.adjust(fisher.smokingtmb[,1],method='BH')

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

fisher.smoking.all <- matrix(0,nrow=dim(mm)[1],ncol=3)
rownames(fisher.smoking.all) <- rownames(mm)
for(i in 1:dim(mm)[1]){
  temp <- fisher.test(rbind(c(sum(mm[i,names.temp1] == 1),
                              sum(mm[i,names.temp1] == 0)),
                            c(sum(mm[i,names.temp2] == 1),
                              sum(mm[i,names.temp2] == 0))))
  fisher.smoking.all[i,1] <- temp$p.value
  fisher.smoking.all[i,3] <- temp$estimate[1]
}
fisher.smoking.all[,2] <- p.adjust(fisher.smoking.all[,1],method='BH')

if(u == 2){
  plot(fisher.smokingtmb[,3],-log10(fisher.smokingtmb[,2]),xlim=c(0,2),ylim=c(0,5),pch=19,cex=2,col='skyblue')
  abline(v=1,h=1,lwd=2)
} else if(u == 3){
  plot(fisher.smokingtmb[,3],-log10(fisher.smokingtmb[,2]),xlim=c(0,2),ylim=c(0,5),pch=19,cex=2,col='red')
  abline(v=1,h=1,lwd=2)
}

##################################################################################################################################
#rna differential expression in TMB high/low & smoking/non-smoking

ttest.smoking.all <- list()
for(j in 1:4){
  if(j == 1){
    names.temp1 <- intersect(intersect(rownames(smoking)[smoking[,3] == 1],names(patient.mutburden)),colnames(rna))
    names.temp1 <- names.temp1[patient.mutburden[names.temp1] < 150]
    names.temp2 <- intersect(intersect(rownames(smoking)[smoking[,3] == 1],names(patient.mutburden)),colnames(rna))
    names.temp2 <- names.temp2[patient.mutburden[names.temp2] >= 150]
  } else if(j == 2){
    names.temp1 <- intersect(intersect(rownames(smoking)[smoking[,3] > 1],names(patient.mutburden)),colnames(rna))
    names.temp1 <- names.temp1[patient.mutburden[names.temp1] < 150]
    names.temp2 <- intersect(intersect(rownames(smoking)[smoking[,3] > 1],names(patient.mutburden)),colnames(rna))
    names.temp2 <- names.temp2[patient.mutburden[names.temp2] >= 150]
  } else if(j == 3){
    names.temp1 <- intersect(intersect(rownames(smoking)[smoking[,3] == 1],names(patient.mutburden)),colnames(rna))
    names.temp1 <- names.temp1[patient.mutburden[names.temp1] < 150]
    names.temp2 <- intersect(intersect(rownames(smoking)[smoking[,3] > 1],names(patient.mutburden)),colnames(rna))
    names.temp2 <- names.temp2[patient.mutburden[names.temp2] < 150]
  } else if(j == 4){
    names.temp1 <- intersect(intersect(rownames(smoking)[smoking[,3] == 1],names(patient.mutburden)),colnames(rna))
    names.temp1 <- names.temp1[patient.mutburden[names.temp1] >= 150]
    names.temp2 <- intersect(intersect(rownames(smoking)[smoking[,3] > 1],names(patient.mutburden)),colnames(rna))
    names.temp2 <- names.temp2[patient.mutburden[names.temp2] >= 150]
  }
  
  names.temp <- intersect(intersect(rownames(smoking)[smoking[,3] > 1],names(patient.mutburden)),colnames(rna))
  ttest.smoking <- matrix(0,nrow=dim(rna)[1],ncol=5)
  rownames(ttest.smoking) <- rownames(rna)
  for(i in 1:dim(rna)[1]){
    temp <- t.test(rna[i,names.temp1],rna[i,names.temp2])
    ttest.smoking[i,1] <- temp$statistic
    ttest.smoking[i,2] <- temp$estimate[1]
    ttest.smoking[i,3] <- temp$estimate[2]
    ttest.smoking[i,4] <- temp$p.value
  }
  ttest.smoking[,5] <- p.adjust(ttest.smoking[,4],method='BH')
  print(j)
  ttest.smoking.all[[j]] <- ttest.smoking
}

plot(ttest.smoking.all[[3]][,1],ttest.smoking.all[[4]][,1])


##################################################################################################################################

#make heatmap for the mutation matrix
mm.pheat <- matrix(0,nrow=150,ncol=dim(mm.sig)[2])
rownames(mm.pheat) <- rownames(mm.sig)
mm.pheat[mm.sig.mis > 0] <- 1
mm.pheat[cna < -2] <- 2
mm.pheat[mm.sig > 0] <- 2
mm.pheat.path <- matrix(0,nrow=17,ncol=dim(mm.sig)[2])
rownames(mm.pheat.path) <- path.names
for(i in 1:17){
  if(sum(genes.sig[,1] == path.names[i]) > 1){
    mm.pheat.path[i,colSums(mm.pheat[genes.sig[,1] == path.names[i],] > 1) > 0] <- 2
    mm.pheat.path[i,colSums(mm.pheat[genes.sig[,1] == path.names[i],] == 1) > 0] <- 1
  } else {
    mm.pheat.path[i,mm.pheat[genes.sig[,1] == path.names[i],] > 1] <- 2
    mm.pheat.path[i,mm.pheat[genes.sig[,1] == path.names[i],] == 1] <- 1
  }
}
pheatmap(mm.pheat.path[order(rowSums(mm.pheat.path > 0),decreasing=TRUE),order(colSums(mm.pheat.path > 0),decreasing=TRUE)],show_colnames = FALSE,cluster_rows = FALSE,cluster_cols = FALSE)
pheatmap(mm.pheat.path[order(rowSums(mm.pheat.path > 0),decreasing=TRUE),c(order(colSums(mm.pheat.path[,1:487] > 0),decreasing=TRUE),487+order(colSums(mm.pheat.path[,488:974] > 0),decreasing=TRUE))],show_colnames = FALSE,cluster_rows = FALSE,cluster_cols = FALSE)

pheatmap(mm.pheat[order(rowSums(mm.pheat > 0),decreasing=TRUE),order(colSums(mm.pheat > 0),decreasing=TRUE)],show_colnames = FALSE,cluster_rows = FALSE,cluster_cols = FALSE)



##################################################################################################################################
#this analysis loads in data that contains the clonality mutations
#Data necessary:
u  = 2
#2 == LUAD, 3 == LUSC

#compile mutation and VAF data
genes.sig <- fread('/Users/michaelsharpnack/Desktop/Kai/DNA%20repair%20related%20genes%20and%20DNA%20polymerases%20genes.csv')
genes.symbol <- genes.sig[[2]]
genes.sig <- as.matrix(genes.sig)
rownames(genes.sig) <- genes.sig[,2]
luad.maf <- fread('/Users/michaelsharpnack/Downloads/tcga_pancancer_dcc_mafs_082115/mafs/tcga_luad_from_dcc.maf',select=c(1,9,11,12,13,16,80,81),skip=195093)
luad.maf2 <- fread('/Users/michaelsharpnack/Downloads/tcga_pancancer_dcc_mafs_082115/mafs/tcga_luad_from_dcc.maf',select=c(1,9,11,12,13,16,80,81),nrows=195091)
luad.maf.3 <- rbind(as.matrix(luad.maf2),as.matrix(luad.maf))
mut.maf <- luad.maf.3
rm(luad.maf2,luad.maf.3,luad.maf)
#list of patients present in both purity and mutation data
tumor.purity <- tumor.purity[intersect(colnames(cna),rownames(tumor.purity)),]
#calculate VAF for each mutation in mut.maf
mut.maf <- cbind(mut.maf,as.numeric(mut.maf[,7])/(as.numeric(mut.maf[,7])+as.numeric(mut.maf[,8])))
mut.maf <- cbind(mut.maf,vector(mode='numeric',dim(mut.maf)[1]))
for(i in 1:dim(mut.maf)[1]){
  if(length(which(rownames(tumor.purity) == mut.maf[i,6])) > 0){
    mut.maf[i,10] <- as.numeric(mut.maf[i,7])/(tumor.purity[mut.maf[i,6],5]*(as.numeric(mut.maf[i,8])-as.numeric(mut.maf[i,8])*(1-tumor.purity[mut.maf[i,6],5])))
  }
  if(i %% 10000 == 0){print(i)}}
genes.symbol.int <- intersect(genes.symbol,unique(mut.maf[,1]))
genes.sig <- genes.sig[genes.symbol.int,]
path.names <- unique(genes.sig[,1])
functional.muts <- union(union(union(union(which(mut.maf[,2] == 'Frame_Shift_Del'),which(mut.maf[,2] == 'Nonsense_Mutation')),
                                     which(mut.maf[,2] == 'Frame_Shift_Ins')),
                               which(mut.maf[,2] == 'In_Frame_Del')),
                         which(mut.maf[,2] == 'In_Frame_Ins'))

mm <- matrix(0,nrow=length(unique(mut.maf[,1])),ncol=length(unique(mut.maf[,6])))
genes.all <- unique(mut.maf[,1])
rownames(mm) <- genes.all
colnames(mm) <- unique(mut.maf[,6])
for(i in 1:length(genes.all)){
  mm[i,unique(mut.maf[,6][which(mut.maf[,1] == genes.all[i])])] <- 1
  if(i %% 1000 == 0){print(i)}
}
colnames(mm) <- substr(colnames(mm),1,15)
class(mm) <- 'numeric'

mm.sig <- matrix(0,nrow=length(genes.symbol.int),ncol=length(unique(mut.maf[,6])))
rownames(mm.sig) <- genes.symbol.int
colnames(mm.sig) <- unique(mut.maf[,6])
for(i in 1:length(genes.symbol.int)){
  mm.sig[i,unique(mut.maf[,6][intersect(which(mut.maf[,1] == genes.symbol.int[i]),functional.muts)])] <- 1
}
colnames(mm.sig) <- substr(colnames(mm.sig),1,15)

#load in copy number alterations from cbioportal
cna.files <- dir('/Users/michaelsharpnack/Desktop/Kai/gistic_cna')
cna <- fread(paste('/Users/michaelsharpnack/Desktop/Kai/gistic_cna/',cna.files[u],sep=''))
#cna <- rbind(cna,fread('/Users/michaelsharpnack/Desktop/Kai/luad_tcga_luad_tcga_gistic.2.txt'))
cna.genes <- cna$GENE_ID
cna <- as.matrix(cna[,3:dim(cna)[2]])
rownames(cna) <- cna.genes
cna <- cna[,colSums(is.nan(cna)) == 0]
mm.sig <- mm.sig[,intersect(colnames(mm.sig),colnames(cna))]
cna <- cna[,intersect(colnames(mm.sig),colnames(cna))]
mm.sig <- mm.sig[intersect(rownames(mm.sig),rownames(cna)),]
cna <- cna[intersect(rownames(mm.sig),rownames(cna)),]
genes.sig <- genes.sig[intersect(rownames(mm.sig),rownames(cna)),]

#calculate each tumor's mutation burden
mut.maf[,6] <- substr(mut.maf[,6],1,15)
patient.mutburden <- vector(mode='numeric',length(unique(mut.maf[,6])))
for(j in 1:length(unique(mut.maf[,6]))){
  patient.mutburden[j] <- length(which(mut.maf[,6] == unique(mut.maf[,6])[j]))
}
names(patient.mutburden) <- unique(mut.maf[,6])
patient.mutburden <- patient.mutburden[colnames(cna)]
#calculate each tumors clonal mutation burden
clonal.mutburden <- vector(mode='numeric',length(unique(mut.maf[,6])))
for(j in 1:length(unique(mut.maf[,6]))){
  clonal.mutburden[j] <- length(which(mut.maf[,6] == unique(mut.maf[,6])[j] & as.numeric(mut.maf[,10]) > 0.4))
}
names(clonal.mutburden) <- unique(mut.maf[,6])
clonal.mutburden <- clonal.mutburden[colnames(cna)]
#calculate each tumors subclonal mutation burden
subclonal.mutburden <- vector(mode='numeric',length(unique(mut.maf[,6])))
for(j in 1:length(unique(mut.maf[,6]))){
  subclonal.mutburden[j] <- length(which(mut.maf[,6] == unique(mut.maf[,6])[j] & as.numeric(mut.maf[,10]) < 0.4))
}
names(subclonal.mutburden) <- unique(mut.maf[,6])
subclonal.mutburden <- subclonal.mutburden[colnames(cna)]

#calculate the mutation burden p-values for sig genes
mut.burden.sig <- matrix(0,nrow=dim(cna)[1],ncol=13)
rownames(mut.burden.sig) <- rownames(cna)
for(i in 1:dim(cna)[1]){
  samples.mut <- colnames(mm.sig)[mm.sig[i,] == 1]
  samples.del <- colnames(mm.sig)[cna[i,] < -1]
  mut.burden.sig[i,1] <- length(union(samples.mut,samples.del))
  mut.burden.sig[i,2] <- mean(patient.mutburden[union(samples.mut,samples.del)],na.rm=TRUE)
  mut.burden.sig[i,3] <- mean(clonal.mutburden[union(samples.mut,samples.del)],na.rm=TRUE)
  mut.burden.sig[i,4] <- mean(subclonal.mutburden[union(samples.mut,samples.del)],na.rm=TRUE)
  #mut.burden.sig[i,5] <- mean(neoantigens.total.filt[union(samples.mut,samples.del)],na.rm = TRUE)
  mut.burden.sig[i,6] <- try(t.test(patient.mutburden[union(samples.mut,samples.del)],patient.mutburden[setdiff(names(patient.mutburden),union(samples.mut,samples.del))])$p.value,silent = TRUE)
  mut.burden.sig[i,8] <- try(t.test(clonal.mutburden[union(samples.mut,samples.del)],clonal.mutburden[setdiff(names(clonal.mutburden),union(samples.mut,samples.del))])$p.value,silent = TRUE)
  mut.burden.sig[i,10] <- try(t.test(subclonal.mutburden[union(samples.mut,samples.del)],subclonal.mutburden[setdiff(names(subclonal.mutburden),union(samples.mut,samples.del))])$p.value,silent = TRUE)
  #mut.burden.sig[i,12] <- try(t.test(neoantigens.total.filt[union(samples.mut,samples.del)],neoantigens.total.filt[setdiff(names(neoantigens.total.filt),union(samples.mut,samples.del))])$p.value,silent = TRUE)
}
mut.burden.sig[,5][nchar(mut.burden.sig[,5])> 100] <- NA
class(mut.burden.sig) <- 'numeric'
#mut.burden.sig <- mut.burden.all[is.na(mut.burden.all[,3]) == FALSE,]
#mut.burden.sig <- mut.burden.all[mut.burden.all[,1] >= 5,]
mut.burden.sig[,7] <- p.adjust(mut.burden.sig[,6],method='BH')
mut.burden.sig[,9] <- p.adjust(mut.burden.sig[,8],method='BH')
mut.burden.sig[,11] <- p.adjust(mut.burden.sig[,10],method='BH')
mut.burden.sig[,13] <- p.adjust(mut.burden.sig[,12],method='BH')
write.csv(mut.burden.sig,'/Users/michaelsharpnack/Desktop/Kai/luad.mutvaf.csv')

################################################################################################################################################
#VAF code 

mut.persiggene <- matrix(0,nrow=length(unique(mut.maf$Hugo_Symbol)),ncol=length(unique(mut.maf$Tumor_Sample_Barcode)))
colnames(mut.persiggene) <- unique(mut.maf$Tumor_Sample_Barcode)
rownames(mut.persiggene) <- unique(mut.maf$Hugo_Symbol)
mut.vaf <- mut.persiggene
temp1 <- as.numeric(levels(mut.maf$V10))[as.numeric(mut.maf$V10)]
for(i in 1:length(unique(mut.maf$Hugo_Symbol))){
  temp <- intersect(which(mut.maf$Hugo_Symbol == unique(mut.maf$Hugo_Symbol)[i]),functional.muts)
  if(length(temp) > 0){
    for(j in 1:length(temp)){
      #mut.persiggene[i,mut.maf$Tumor_Sample_Barcode[temp[j]]] <- sum(temp1[which(mut.maf$Tumor_Sample_Barcode == mut.maf$Tumor_Sample_Barcode[temp[j]])] <= temp1[temp[j]])    
      mut.vaf[i,mut.maf$Tumor_Sample_Barcode[temp[j]]] <- temp1[temp[j]]
    }
  }
  if(i %% 1000 == 0){print(i)}
}
mut.persiggene <- mut.persiggene[,patients.temp]
mut.persiggene[mut.persiggene == 0] <- NA
mut.vaf <- mut.vaf[,patients.temp]
mut.persiggene.vaf <- mut.persiggene/mut.vaf
mut.persiggene.vaf[is.nan(mut.persiggene.vaf)] <- 0
mut.vaf[mut.vaf == 0] <- NA
#calculate each tumor's average & standard deviation of vaf
mut.maf$Tumor_Sample_Barcode <- substr(mut.maf$Tumor_Sample_Barcode,1,15)
patient.vaf <- matrix(0,nrow=length(unique(mut.maf$Tumor_Sample_Barcode)),ncol=2)
for(j in 1:length(unique(mut.maf$Tumor_Sample_Barcode))){
  patient.vaf[j,1] <- mean(as.numeric(levels(mut.maf$V10))[mut.maf$V10][which(mut.maf$Tumor_Sample_Barcode == unique(mut.maf$Tumor_Sample_Barcode)[j])])
  patient.vaf[j,2] <- sd(as.numeric(levels(mut.maf$V10))[mut.maf$V10][which(mut.maf$Tumor_Sample_Barcode == unique(mut.maf$Tumor_Sample_Barcode)[j])])
}
rownames(patient.vaf) <- unique(mut.maf$Tumor_Sample_Barcode)
patient.vaf <- patient.vaf[patients.temp,]

