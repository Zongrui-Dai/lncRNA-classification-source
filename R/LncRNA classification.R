library(dplyr)
library(mlbench)
library(caret)
library(glmnet)
library(pROC)
library(gmodels)
library(h2o)
library(reshape2)
library(ggstatsplot)
library(e1071)
library(elmNNRcpp)
setwd('D:/test/lncRNA/')
lncRNA<-read.csv('lncRNA_training.csv')
cRNA<-read.csv('codingRNA_training.csv')

Coding <- cRNA[,c(-1,-13)]
Noncoding <- lncRNA[,c(-1,-13)]
dataset<-rbind(Coding,Noncoding)

d<-max_min_scale(dataset,dataframe = T)
#####################################################
y = as.matrix(dataset$label)
x = as.matrix(dataset[,-length(dataset)])
fit1<-cv.glmnet(x=x,y=y,alpha = 1,family = 'binomial',type.measure="auc",nfolds = 10)
plot(fit1)

lasso<-data.frame(fit1$glmnet.fit$df,fit1$cvm,fit1$lambda)
#####################################################
id<-c();remain<-c()
for (i in seq(0.1,1,0.1)){
  remove_trans <- findCorrelation(cor(dataset[,-length(dataset)]),cutoff=i)
  remain<-c(remain,11-length(remove_trans))
  id<-c(id,i)
}
pearson<-data.frame(id,remain)

#####################################################
pca_auc<-c();lasso_auc<-c();pearson_auc<-c()
for(i in 11){
  cvlist <- CVgroup(k = 10,datasize = nrow(dataset),seed = 123)
  t<-cvlist[[i-1]]
  
  d.pca<-princomp(dataset[,-length(dataset)],cor=T,scores=T) 
  pca_dataset<-d.pca$scores[,1:i]
  pca<-svm(x=as.matrix(pca_dataset[-t,]),y=as.matrix(dataset[-t,'label']),kernel = 'radial')
  pca_pre<-predict(pca,as.matrix(pca_dataset[t,]))
  r<-roc(dataset[t,'label'],pca_pre,auc=T,thres=TRUE)
  pca_auc<-c(pca_auc,r$auc)
  
  l1<-lasso[lasso$fit1.glmnet.fit.df == i,]
  lambda<-l1[l1$fit1.cvm == max(l1$fit1.cvm),3]
  param<-coef(fit1,s=lambda)
  param<-as.data.frame(as.matrix(param))
  param$feature<-rownames(param)
  param_e2<-param[param$'1'!= 0,]
  lasso_feature<-rownames(param_e2[-1])[-1]
  las<-dataset[,lasso_feature]
  lass<-svm(x=as.matrix(las[-t,]),y=as.matrix(dataset[-t,'label']),kernel = 'radial')
  lass_pre<-predict(lass,as.matrix(las[t,]))
  r<-roc(dataset[t,'label'],lass_pre,auc=T,thres=TRUE)
  lasso_auc<-c(lasso_auc,r$auc)
  
  p<-pearson[pearson$remain==i,]
  p<-p[1,1]
  remove_trans <- findCorrelation(cor(dataset[,-length(dataset)]),cutoff=p)
  train_trans <- dataset[,-remove_trans]
  pear<-svm(x=as.matrix(train_trans[-t,]),y=as.matrix(dataset[-t,'label']),kernel = 'radial')
  pear_pre<-predict(pear,as.matrix(train_trans[t,]))
  r<-roc(dataset[t,'label'],pear_pre,auc=T,thres=TRUE)
  pearson_auc<-c(pearson_auc,r$auc)
}
##########################################################
l1<-lasso[lasso$fit1.glmnet.fit.df == 8,]
lambda<-l1[l1$fit1.cvm == max(l1$fit1.cvm),3]
param<-coef(fit1,s=lambda)
param<-as.data.frame(as.matrix(param))
param$feature<-rownames(param)
param_e2<-param[param$'1'!= 0,]
lasso_feature<-rownames(param_e2[-1])[-1]

a2<-dataset[,lasso_feature]
datasize <- nrow(a2)
k=10
cvlist <- CVgroup(k = k,datasize = datasize,seed = 123)
h2o.init()
b1<-cbind(dataset$label,a2)
colnames(b1)[1]<-'label'

sen<-c();spe<-c();ac<-c();auc<-c()
for(i in 1:k){
  train <- as.data.frame(b1[-cvlist[[i]],])  
  test <- as.data.frame(b1[cvlist[[i]],])
  train1<-as.h2o(train)
  rf<-h2o.randomForest(y=1,training_frame = train1,nfolds = k,ntrees = 50)
  p1<-predict(rf,as.h2o(test))
  p<-as.data.frame(p1)
  r<-roc(test$label,p$predict,auc=T,thres=TRUE)
  e<-cbind(r$thresholds,r$sensitivities+r$specificities)
  best<-subset(e,e[,2]==max(e[,2]))[,1]
  best<-as.numeric(best)
  p[p$predict<best,]<-0
  p[p$predict>=best,]<-1
  a<-CrossTable(test$label,p$predict,prop.c = F,prop.t = F,prop.chisq = F)
  sen<-c(sen,a$t[1,1]/(a$t[1,1]+a$t[1,2]))
  spe<-c(spe,a$t[2,2]/(a$t[2,2]+a$t[2,1]))
  ac<-c(ac,(a$t[1]+a$t[4])/sum(a$t))
  p11<-predict(rf,as.h2o(train))
  auc<-c(auc,r$auc)
  
  rf<-h2o.gbm(y=1,training_frame = train1,nfolds = k,ntrees = 50)
  p2<-predict(rf,as.h2o(test))
  p<-as.data.frame(p2)
  r<-roc(test$label,p$predict,auc=T,thres=TRUE)
  e<-cbind(r$thresholds,r$sensitivities+r$specificities)
  best<-subset(e,e[,2]==max(e[,2]))[,1]
  best<-as.numeric(best)
  p[p$predict<best,]<-0
  p[p$predict>=best,]<-1
  a<-CrossTable(test$label,p$predict,prop.c = F,prop.t = F,prop.chisq = F)
  sen<-c(sen,a$t[1,1]/(a$t[1,1]+a$t[1,2]))
  spe<-c(spe,a$t[2,2]/(a$t[2,2]+a$t[2,1]))
  ac<-c(ac,(a$t[1]+a$t[4])/sum(a$t))
  p12<-predict(rf,as.h2o(train))
  auc<-c(auc,r$auc)
  
  rf<-elm_train(x = as.matrix(train[,-1]),y=as.matrix(train[,1]),nhid = 50, actfun = 'relu')
  p3<-elm_predict(rf,as.matrix(test[,-1]))
  p<-as.data.frame(p3)
  r<-roc(test$label,p$V1,auc=T,thres=TRUE)
  e<-cbind(r$thresholds,r$sensitivities+r$specificities)
  best<-subset(e,e[,2]==max(e[,2]))[,1]
  best<-as.numeric(best)
  p[p$V1<best,]<-0
  p[p$V1>=best,]<-1
  a<-CrossTable(test$label,p$V1,prop.c = F,prop.t = F,prop.chisq = F)
  sen<-c(sen,a$t[1,1]/(a$t[1,1]+a$t[1,2]))
  spe<-c(spe,a$t[2,2]/(a$t[2,2]+a$t[2,1]))
  ac<-c(ac,(a$t[1]+a$t[4])/sum(a$t))
  p13<-elm_predict(rf,as.matrix(train[,-1]))
  auc<-c(auc,r$auc)
  
  rf<-h2o.glm(y=1,training_frame = train1,nfolds = k,alpha = 1)
  p4<-predict(rf,as.h2o(test))
  p<-as.data.frame(p4)
  r<-roc(test$label,p$predict,auc=T,thres=TRUE)
  e<-cbind(r$thresholds,r$sensitivities+r$specificities)
  best<-subset(e,e[,2]==max(e[,2]))[,1]
  best<-as.numeric(best)
  p[p$predict<best,]<-0
  p[p$predict>=best,]<-1
  a<-CrossTable(test$label,p$predict,prop.c = F,prop.t = F,prop.chisq = F)
  sen<-c(sen,a$t[1,1]/(a$t[1,1]+a$t[1,2]))
  spe<-c(spe,a$t[2,2]/(a$t[2,2]+a$t[2,1]))
  ac<-c(ac,(a$t[1]+a$t[4])/sum(a$t))
  p14<-predict(rf,as.h2o(train))
  auc<-c(auc,r$auc)
  
  p11<-as.data.frame(p11);p12<-as.data.frame(p12);p13<-as.data.frame(p13);p14<-as.data.frame(p14)
  m<-data.frame(train$label,p11,p12,p13,p14);m<-as.h2o(m)
  n<-data.frame(test$label,as.data.frame(p1),as.data.frame(p2),as.data.frame(p3),as.data.frame(p4));n<-as.h2o(n)
  sec<-h2o.randomForest(y=1,training_frame = m,nfolds = k,,ntrees = 100)
  p5<-predict(sec,as.h2o(n))
  p<-as.data.frame(p5)
  r<-roc(test$label,p$predict,auc=T,thres=TRUE)
  e<-cbind(r$thresholds,r$sensitivities+r$specificities)
  best<-subset(e,e[,2]==max(e[,2]))[,1]
  best<-as.numeric(best)
  p[p$predict<best,]<-0
  p[p$predict>=best,]<-1
  a<-CrossTable(test$label,p$predict,prop.c = F,prop.t = F,prop.chisq = F)
  sen<-c(sen,a$t[1,1]/(a$t[1,1]+a$t[1,2]))
  spe<-c(spe,a$t[2,2]/(a$t[2,2]+a$t[2,1]))
  ac<-c(ac,(a$t[1]+a$t[4])/sum(a$t))
  auc<-c(auc,r$auc)
}
library('cowplot')
type<-c(rep('rf',k),rep('gbm',k),rep('elm',k),rep('glm',k),rep('mix',k))
g<-c(1:k,1:k,1:k,1:k,1:k)
data<-data.frame(type,sen,spe,ac,auc,g)
p1<-ggplot(data,aes(x=g,y=sen,group=type))+geom_line(aes(colour=type))+geom_point(aes(colour=type))+ggtitle('Sensitive ability')+geom_text(aes(label=round(sen,3)), size=4)
p2<-ggplot(data,aes(x=g,y=spe,group=type))+geom_line(aes(colour=type))+geom_point(aes(colour=type))+ggtitle('Specific ability')+geom_text(aes(label=round(spe,3)), size=4)
p3<-ggplot(data,aes(x=g,y=ac,group=type))+geom_line(aes(colour=type))+geom_point(aes(colour=type))+ggtitle('Accuarcy ability')+geom_text(aes(label=round(ac,3)), size=4)
p4<-ggplot(data,aes(x=g,y=auc,group=type))+geom_line(aes(colour=type))+geom_point(aes(colour=type))+ggtitle('AUC')+geom_text(aes(label=round(sen,3)), size=4)
plot_grid(p1,p2,p3,p4)

p1<-ggplot(data,aes(x=type,y=sen,group=type))+geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+ggtitle('Sensitive ability')
p2<-ggplot(data,aes(x=type,y=spe,group=type))+geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+ggtitle('Specific ability')
p3<-ggplot(data,aes(x=type,y=ac,group=type))+geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+ggtitle('Accuarcy ability')
p4<-ggplot(data,aes(x=type,y=auc,group=type))+geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+ggtitle('AUC')
plot_grid(p1,p2,p3,p4)


####################################################
##Calulate the max_min scale
max_min_scale<-function(data,dataframe=F){
  if(dataframe==F){
    data<-(data-min(data))/(max(data)-min(data))
  }else{
    for(i in 1:length(data)){
      scale<-(data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
      data[,i]<-scale
    }
    return(data)
  }
}
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize] 
  temp <- sample(n,datasize)   
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])
}

p1<-ggplot(data,aes(x=type,y=auc,fill=type))+ 
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ 
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+ggtitle('AUC')
p2<-ggplot(data,aes(x=type,y=spe,fill=type))+ 
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ 
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+ggtitle('Specifity')
p3<-ggplot(data,aes(x=type,y=sen,fill=type))+ 
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ 
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+ggtitle('Sensitivity')
p4<-ggplot(data,aes(x=type,y=ac,fill=type))+ 
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ 
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+ggtitle('Accuracy')
plot_grid(p1,p2,p3,p4)

