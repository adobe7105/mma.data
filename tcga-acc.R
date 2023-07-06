library(limma)   
library(dplyr)
expFile="TCGA-BLCA-FPKM.symbol.txt"       
geneFile="gene.txt"     
setwd("D:\\2.mma\\提取与mma代谢相关的基因\\TCGA-BLCA")
###########挑选基因###########
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
gene<-read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]
geneExp=log2(geneExp+1)
Tab=rbind(ID=colnames(geneExp),geneExp)
write.table(Tab, file="TCGA-BLCA-mmaExp.txt", sep="\t", quote=F, col.names=F)
##读取并挑选我们想要的基因的表达量

##########id处理#############
#############################
#############################
#############################
lncFile="TCGA-BLCA-mmaExp.txt"  
cliFile="TCGA-BLCA-clinical.csv"
rt.m=read.table(lncFile, header=T, sep="\t", check.names=F,row.names = 1)
rt.m=t(as.matrix(rt.m))
rt.m=avereps(rt.m)
group=sapply(strsplit(rownames(rt.m),"\\."), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rownames(rt.m)=gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*", "\\1\\.\\2\\.\\3", rownames(rt.m))
rt.tumor= rt.m[group==0,]
rt.normal=rt.m[group==1,]

######筛选肿瘤样本#########
###########################
###########################
rt.clinical <- read.csv("D:/2.mma/提取与mma代谢相关的基因/TCGA-BLCA/TCGA-BLCA-clinical.csv", row.names=1)
colnames(rt.clinical)#挑选一下要哪些
aa=c("disease","days_to_death","year_of_death","submitter_id",
     "gender" ,"age_at_diagnosis","age_at_index","vital_status",
     "days_to_last_follow_up","tumor_stage"
)
#在这里，对于存活的样本，我们采用最后随访时间作为失访时间，对于死亡样本以死亡时间作为随访时间，以确证年龄作为
#样本年龄，我强烈建议别全部按照我的来，因为我是直接抄的整理好的数据，不同肿瘤的i细节不一样，你应该用tcgabiolinks下载的为准。
rt.clinical=rt.clinical[,aa]
rt.clinical.dead=rt.clinical %>% filter(vital_status=="Dead")
rt.clinical.clean.dead=cbind(id=rt.clinical.dead$submitter_id,
                             os=rt.clinical.dead$vital_status,
                             time=(rt.clinical.dead$days_to_death)/365,
                             age=(rt.clinical.dead$age_at_diagnosis)/365
                             #,stage=rt.clinical.dead$tumor_stage
                             #,gender=rt.clinical.dead$gender
)
rt.clinical.alive=rt.clinical %>% filter(vital_status=="Alive")
rt.clinical.clean.alive=cbind(id=rt.clinical.alive$submitter_id,
                              os=rt.clinical.alive$vital_status,
                              time=(rt.clinical.alive$days_to_last_follow_up)/365,
                              age=(rt.clinical.alive$age_at_diagnosis)/365
                              #,stage=rt.clinical.alive$tumor_stage
                              #,gender=rt.clinical.alive$gender
)


rt.clinical.clean=as.data.frame(rbind(rt.clinical.clean.alive,rt.clinical.clean.dead))
rt.clinical.clean$id=gsub("-",".",as.vector(rt.clinical.clean[,1]))
write.table(rt.clinical.clean, file="TCGA-BLCA-clinical-clean.txt", sep="\t", quote=F, col.names=T)

######合并临床数据和表达数据####
cl=read.table("TCGA-BLCA-clinical-clean.txt", header=T, sep="\t", check.names=F,row.names = 2)
samesample=intersect(rownames(cl), rownames(rt.tumor))
exp=rt.tumor[samesample,]
cl=cl[samesample,]
all1=cbind(cl[,-1],exp)

samesample2=intersect(rownames(cl), rownames(rt.normal))
exp2=rt.normal[samesample2,]
cl=cl[samesample2,]
all2=cbind(cl[,-1],exp2)
###########################################################

old.tumor=all1%>%filter(age>=60)
young.tumor=all1%>% filter(age<=30)

write.table(all1, file="TCGA-BLCA-merge.tumor.txt", sep="\t", quote=F, col.names=T)
write.table(all2, file="TCGA-BLCA-merge.normal.txt", sep="\t", quote=F, col.names=T)
write.table(old.tumor, file="TCGA-BLCA-merge.tumor.old.tumer.txt", sep="\t", quote=F, col.names=T)
write.table(young.tumor, file="TCGA-BLCA-merge.tumor.young.tumor.txt", sep="\t", quote=F, col.names=T)



##########################################
##########################################
##########################################
##########################################
##########################################



######################################################################

#单因素cox
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
library(DynNom)
all1$os=as.numeric(as.factor(all1$os))-1 
all1=all1%>%filter(time>0)
all1=all1%>%filter(age>0)

sigGenes<-c("os","time")
outUniTab=data.frame()
for(i in colnames(all1[,4:ncol(all1)])){
  if(sd(all1[,i])>0.1){
    
    cox <- coxph(Surv(time,os) ~ all1[,i], data = all1)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    
    if(coxP<0.05){
      sigGenes=c(sigGenes,i)
      outUniTab=rbind(outUniTab,
                      cbind(id=i,
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            HR.95L=coxSummary$conf.int[,"lower .95"],
                            HR.95H=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
  }	
}
#提取基因
uniSigExp=all1[,sigGenes]
uniSigExp=na.omit(uniSigExp)
uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)

#lasso
set.seed(1234)
x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
y=data.matrix(Surv(uniSigExp$time,uniSigExp$os))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
coef <- coef(cvfit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoSigExp=uniSigExp[,c("os", "time", lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
#提取基因lassosigexp

multiCox <- coxph(Surv(time,os) ~ ., data = lassoSigExp)
multiCox=step(multiCox, direction = "both")
multiCoxSum=summary(multiCox)

#森林图准备
outMultiTab=data.frame()
outMultiTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
HR=as.numeric(outMultiTab[,2])
HR=as.data.frame(HR)
gene=rownames(outMultiTab)
gene=as.data.frame(gene)
outMultiTab=cbind(gene=gene,HR=HR)
outMultiTab=outMultiTab%>%arrange(HR)%>%mutate(gene=factor(gene,levels=gene))
#############################################
riskScore=predict(multiCox,type="risk",newdata=all1) 
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("OS","Censor",coxGene)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
testRiskOut=cbind(all1,riskScore,risk)
fit1 <- survfit(Surv(time ,os) ~risk, data = testRiskOut)
##第一张图，与甲基丙二酸代谢相关的基因所制作的预测模型，lasso+cox确定危险因素

pdf("TCGA-BLCA-coxmodel-lasso.pdf",width=5,heigh=5)
p1<-ggsurvplot(fit1,                     # survfit object with calculated statistics.
               title="In all sample",
               pval = TRUE,             # show p-value of log-rank test.
               conf.int = TRUE,         # show confidence intervals for 
               # point estimaes of survival curves.
               conf.int.style = "step",  # customize style of confidence intervals
               xlab = "Time in years",   # customize X axis label.
               ggtheme = theme_light(), # customize plot and risk table with a theme.
               #risk.table = "abs_pct",  # absolute number and percentage at risk.
               risk.table.y.text.col = T,# colour risk table text annotations.
               risk.table.y.text = FALSE,# show bars instead of names in text annotations
               surv.median.line = "hv",  # add the median survival pointer.
               palette =  c("#E7B800", "#2E9FDF") # custom color palettes.
)
p1
dev.off()

p7<-ggplot(outMultiTab,aes(x=gene,y=HR))+geom_segment(aes(x=gene,xend=gene,y=0,yend=HR),
                                                      size=1.5,color="#C9CACA",linetype="solid")+
  geom_point(size=5,color="#8e8bff",fill="#fea3a2",shape=21)+theme_light()+
  theme(panel.grid.major.x=element_blank(),
        panel.border=element_blank(),
        axis.ticks.x=element_blank())+
  xlab("")+ylab("HR in cox regression")+coord_flip()
p7

pdf("TCGA-BLCA-BBT.pdf",width=4,heigh=4)
p8<-ggplot(outMultiTab,aes(x=gene,y=HR))+geom_segment(aes(x=gene,xend=gene,y=1,yend=HR), size=ifelse(outMultiTab$gene%in%c("MMAA","MMUT","MUT","MMAB"),1,2), color=ifelse(outMultiTab$gene%in%c("MMAA","MMUT","MUT","MMAB"),"#F06955","#C9CACA"))+
  geom_point(size=ifelse(outMultiTab$gene%in%c("MMAA","MMUT","MUT","MMAB"),8,5),
             color=ifelse(outMultiTab$gene%in%c("MMAA","MMUT","MUT","MMAB"),"#f6b37f","#74C6BE"), fill=ifelse(outMultiTab$gene%in%c("MMAA","MMUT","MUT","MMAB"),"#F06955","#74C6BE"))+
  theme_light()+theme(panel.grid.major.x=element_blank(),
                      panel.border=element_blank(),
                      axis.ticks.x=element_blank())+xlab("")+ylab("HR in cox regression")+
  coord_flip()
p8
dev.off()
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
medianTrainRisk2=median(all1$MMAA)
risk_MMAA=as.vector(ifelse(all1$MMAA>medianTrainRisk2,"high_MMAA","low_MMAA"))
testRiskOut2=cbind(all1,risk_MMAA)
fit2 <- survfit(Surv(time ,os) ~risk_MMAA, data = testRiskOut2)
#######################################################################
old.tumor$os=as.numeric(as.factor(old.tumor$os))-1 
old.tumor=old.tumor%>%filter(time>=0)
old.tumor=old.tumor%>%filter(age>0)
risk_MMAA=as.vector(ifelse(old.tumor$MMAA>medianTrainRisk2,"high_MMAA","low_MMAA"))
testRiskOut3=cbind(old.tumor,risk_MMAA)
fit3 <- survfit(Surv(time ,os) ~risk_MMAA, data = testRiskOut3)


###第二章，采用mut基因的预测，这个其实你可以换成其他几个

p2<-ggsurvplot(
  fit2,                     # survfit object with calculated statistics.
  title="Based on MMAA expressions in all samples",
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in years",   # customize X axis label.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  #risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  surv.median.line = "hv",  # add the median survival pointer.
  palette =  c("#96c37d", "#c497b2") # custom color palettes.
)
pdf("TCGA-BLCA-mmaa in all.pdf",width=5,heigh=5)
p2
dev.off()
###################################

p3<-ggsurvplot(
  fit3,                     # survfit object with calculated statistics.
  title="Based on MMAA expression in old samples",
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in years",   # customize X axis label.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  #risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  surv.median.line = "hv",  # add the median survival pointer.
  palette =  c("#8e8bff", "#fea3a2") # custom color palettes.
)
pdf("TCGA-BLCA-mmaa in old.pdf",width=5,heigh=5)
p3
dev.off()
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
medianTrainRisk33=median(all1$MUT)
risk_MUT=as.vector(ifelse(all1$MUT>medianTrainRisk33,"high_MUT","low_MUT"))
testRiskOut33=cbind(all1,risk_MUT)
fit22 <- survfit(Surv(time ,os) ~risk_MUT, data = testRiskOut33)
#######################################################################
old.tumor$os=as.numeric(as.factor(old.tumor$os))-1 
old.tumor=old.tumor%>%filter(time>=0)
old.tumor=old.tumor%>%filter(age>0)
risk_MUT=as.vector(ifelse(old.tumor$MUT>medianTrainRisk33,"high_MUT","low_MUT"))
testRiskOut33=cbind(old.tumor,risk_MUT)
fit33 <- survfit(Surv(time ,os) ~risk_MUT, data = testRiskOut33)


###第二章，采用mut基因的预测，这个其实你可以换成其他几个

p22<-ggsurvplot(
  fit22,                     # survfit object with calculated statistics.
  title="Based on MUT expressions in all samples",
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in years",   # customize X axis label.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  #risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  surv.median.line = "hv",  # add the median survival pointer.
  palette =  c("#96c37d", "#c497b2") # custom color palettes.
)
pdf("TCGA-BLCA-mut in all.pdf",width=5,heigh=5)
p22
dev.off()
###################################

p33<-ggsurvplot(
  fit33,                     # survfit object with calculated statistics.
  title="Based on MUT expression in old samples",
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in years",   # customize X axis label.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  #risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  surv.median.line = "hv",  # add the median survival pointer.
  palette =  c("#8e8bff", "#fea3a2") # custom color palettes.
)
pdf("TCGA-BLCA-mut in old.pdf",width=5,heigh=5)
p33
dev.off()
splots <- list()
splots[[1]] <- p1
splots[[2]] <- p2
# 将多个图合并一起


