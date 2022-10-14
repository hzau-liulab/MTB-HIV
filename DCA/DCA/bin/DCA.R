args=commandArgs(T)
disease_input=args[1]
HC_input=args[2]
ppi=args[3]
DiseaseState=args[4]
out=args[5]

##usage: Rscript DCA.R ../examples/input/HMI ../examples/input/HC ../data/PPI.txt HMI ../examples/output/

library(plyr)
ifelse(!dir.exists(file.path(args[5])), dir.create(file.path(args[5])), FALSE)
dis_express<- dir(path=args[1],pattern="_expression.txt")
if (length(dis_express)!=0){
  for (i in 1:length(dis_express)){
    file_name<-strsplit(as.character(dis_express[i]),'_expression')[[1]][1]
    expression<-read.table(paste(args[1],dis_express[i],sep='/'),header=T,row.names = "ID_REF",sep="\t")
    
	sample_name<-paste(file_name,'_Infor.txt',sep='')
	name<-read.table(paste(args[1],sample_name,sep='/'),header=F,sep="\t")
    
    ptb_name<-name[grep(DiseaseState,name$V2),]$V1
    ptb_exp<-subset(expression,select=as.character(ptb_name))
	
	pcc<-cor(t(ptb_exp),method = "pearson")
	pairs<-read.table(file=ppi,header=F)
	pairs$id<-row.names(pairs)
	pairs_pcc<-ddply(pairs,.(id),function(x){ if(x$V1 %in% row.names(pcc) & x$V2 %in% row.names(pcc)){pcc[as.character(x$V1),as.character(x$V2)]}else{ "NA" }})
	write.table(file=paste(out,"pairs.txt",sep="/"),pairs,quote=F,row.names=F)
	file1<-paste(file_name,'_disease_pcc.txt',sep='')
	write.table(file=paste(out,file1,sep="/"),pairs_pcc,quote=F,row.names=F)
  }
}

hc_express<- dir(path=args[2],pattern="_expression.txt")
if (length(hc_express)!=0){
  for (i in 1:length(hc_express)){
    file_name<-strsplit(as.character(hc_express[i]),'_expression')[[1]][1]
    expression<-read.table(paste(args[2],hc_express[i],sep='/'),header=T,row.names = "ID_REF",sep="\t")
    
	sample_name<-paste(file_name,'_Infor.txt',sep='')
	name<-read.table(paste(args[2],sample_name,sep='/'),header=F,sep="\t")
    con_name<-name[grep("HC",name$V2),]$V1
    con_exp<-subset(expression,select=as.character(con_name))

	pcc<-cor(t(con_exp),method = "pearson")
	pairs<-read.table(file=ppi,header=F)
	pairs$id<-row.names(pairs)
	pairs_pcc<-ddply(pairs,.(id),function(x){ if(x$V1 %in% row.names(pcc) & x$V2 %in% row.names(pcc)){pcc[as.character(x$V1),as.character(x$V2)]}else{ "NA" }})
	write.table(file=paste(out,"pairs.txt",sep="/"),pairs,quote=F,row.names=F)
	file2<-paste(file_name,'_HC_pcc.txt',sep='')
	write.table(file=paste(out,file2,sep="/"),pairs_pcc,quote=F,row.names=F)
  }
}


###
dis_pcc<- dir(path=args[5],pattern="_disease_pcc.txt")
dis_data<-read.table(paste(args[5],dis_pcc[1],sep='/'),header = T,dec = '.')

if (length(dis_pcc)!=0){
  for (i in 2:length(dis_pcc)){
    new.data = read.table(paste(args[5],dis_pcc[i],sep='/'), header=TRUE, dec = ".")
    dis_data = merge(dis_data,new.data,by='id')
  }
}
dis_data[is.na(dis_data)]=0
dis_data=as.data.frame(lapply(dis_data,as.numeric))
dis_data$mean_dis_fc<- rowSums(dis_data[,2:5])/4
dis_data <- dis_data[,c('id','mean_dis_fc')]

hc_pcc<- dir(path=args[5],pattern="_HC_pcc.txt")
hc_data<-read.table(paste(args[5],hc_pcc[1],sep='/'),header = T,dec = '.')
new_hc_data<-read.table(paste(args[5],hc_pcc[2],sep='/'), header=TRUE, dec = ".")
hc_data <- merge(hc_data,new_hc_data,by='id')
hc_data[is.na(hc_data)]=0
hc_data=as.data.frame(lapply(hc_data,as.numeric))
hc_data$mean_hc_fc<- rowSums(hc_data[,2:3])/2
hc_data <- hc_data[,c('id','mean_hc_fc')]

logfc_data<-merge(dis_data,hc_data,by='id')
logfc_data$abs_pcc_cha <- abs(logfc_data$mean_dis_fc-logfc_data$mean_hc_fc)
Thre<-round(round(mean(logfc_data$abs_pcc_cha),3)+3*round(sd(logfc_data$abs_pcc_cha),3),2)
step1 <- subset(logfc_data,logfc_data$abs_pcc_cha >= Thre,select= c(1:4))


step2 <- subset(step1,(abs(step1$mean_dis_fc)>0.5 & abs(step1$mean_hc_fc)<0.5) | (abs(step1$mean_dis_fc)<0.5 & abs(step1$mean_hc_fc)>0.5),select= c(1:4))

step3 <- subset(step2,(step2$mean_dis_fc>0 & step2$mean_hc_fc<0) | (step2$mean_dis_fc<0 & step2$mean_hc_fc>0),select= c(1:4))

pair_file<-read.table(paste(args[5],"pairs.txt",sep='/'),header=T,sep = ' ',col.names = c("pro1","pro2","id"))
pair_file <- pair_file[,c('id','pro1','pro2')]

gene_pair<-merge(step3,pair_file,by='id')
gene_pair <- gene_pair[,c("pro1",'pro2')]
write.table(file=paste(args[5],"gene_pair.txt",sep='/'),gene_pair,quote=F,row.names=F)



