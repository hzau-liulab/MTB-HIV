args=commandArgs(T)
expression=args[1]
sampleInfor=args[2]
DiseaseState=args[3]
out=args[4]

library(limma)
data<-read.table(file=expression,header=T,row.names = "ID_REF")

name<-read.table(file=sampleInfor,header=F,sep="\t")
ptb_name<-name[grep(DiseaseState,name$V2),]$V1        
con_name<-name[grep("HC",name$V2),]$V1

all<-c(ptb_name,con_name)
exp_select<-subset(data,select=all)

samps <- factor(c(rep("ptb",length(ptb_name)),rep("control",length(con_name))))

design <- model.matrix(~0+samps)
rownames(design)=colnames(exp_select)
colnames(design) = c("ptb","control")
cont.matrix<-makeContrasts(control-ptb,levels=design)    
fit <- lmFit(exp_select,design)
fit2 <- contrasts.fit(fit, cont.matrix) 
fit2 <- eBayes(fit2)                             
tempOutltb = topTable(fit2, coef=1, n=Inf)
#tabltb <- topTable(fit2, adjust = "BH", confin =T , number = 1000, sort.by="logFC")

nrDEG_ltb = na.omit(tempOutltb)
nrDEG_ltb$logFC<-abs(nrDEG_ltb$logFC)
Thre<-round(round(mean(nrDEG_ltb$logFC),3)+3*round(sd(nrDEG_ltb$logFC),3),2)
sub_ltb <- subset(nrDEG_ltb,nrDEG_ltb$logFC >= Thre,select= c(1:6))
deg <- subset(sub_ltb,sub_ltb$adj.P.Val < 0.05,select= c(1:6) )
write.table(file=paste(out,"deg.txt",sep = ""),deg,quote = F,row.names = T,sep="\t")