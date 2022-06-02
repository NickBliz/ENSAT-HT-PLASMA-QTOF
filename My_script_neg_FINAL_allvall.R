setwd("D:/Office 2010 Toolkit and EZ-Activator v2.2.3/metabolomics hypertension/QTOF_EQP_ENSAT-HT/ALL_NEG/BLANKS")
EQP=read.table("EQP.tsv", sep="\t", header = T)
EQP_names <- read_csv("EQP_names.csv")
EQP_names=t(EQP_names)
EQP_names=EQP_names[-(1:2),]
colnames(EQP_names)=colnames(EQP)
EQP_b1d1=EQP[,which(EQP_names[10,]=="run1" & (EQP_names[9,]==1))]
EQP_b1d2=EQP[,which(EQP_names[10,]=="run1" & (EQP_names[9,]==2))]
EQP_b2d1=EQP[,which(EQP_names[10,]=="run3" & (EQP_names[9,]==1))]
EQP_b2d2=EQP[,which(EQP_names[10,]=="run3" & (EQP_names[9,]==2))]
EQP_b3d1=EQP[,which(EQP_names[10,]=="run4" & (EQP_names[9,]==1))]
EQP_b3d2=EQP[,which(EQP_names[10,]=="run4" & (EQP_names[9,]==2))]
EQP_b4d1=EQP[,which(EQP_names[10,]=="run2" & (EQP_names[9,]==1))]
EQP_b4d2=EQP[,which(EQP_names[10,]=="run2" & (EQP_names[9,]==2))]
EQP_b5d1=EQP[,which(EQP_names[10,]=="run5" & (EQP_names[9,]==1))]
EQP_b5d2=EQP[,which(EQP_names[10,]=="run5" & (EQP_names[9,]==2))]
EQP_b6d1=EQP[,which(EQP_names[10,]=="run6" & (EQP_names[9,]==1))]
EQP_b6d2=EQP[,which(EQP_names[10,]=="run6" & (EQP_names[9,]==2))]

EQP_b1d1=cbind(EQP[,(1:9)], EQP_b1d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b1d2=cbind(EQP[,(1:9)], EQP_b1d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b2d1=cbind(EQP[,(1:9)], EQP_b2d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b2d2=cbind(EQP[,(1:9)], EQP_b2d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b3d1=cbind(EQP[,(1:9)], EQP_b3d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b3d2=cbind(EQP[,(1:9)], EQP_b3d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b4d1=cbind(EQP[,(1:9)], EQP_b4d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b4d1=EQP_b4d1[,-29]
EQP_b4d2=cbind(EQP[,(1:9)], EQP_b4d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
#REMOVE ALL SAMPLES NOT RERUN, INCLUDING QCs AT BEGINNING!
EQP_b4d2=EQP_b4d2[,-(32:35)]
EQP_b5d1=cbind(EQP[,(1:9)], EQP_b5d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b5d2=cbind(EQP[,(1:9)], EQP_b5d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b6d1=cbind(EQP[,(1:9)], EQP_b6d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b6d2=cbind(EQP[,(1:9)], EQP_b6d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])

bl=read.table("EQP_blanks.tsv", sep="\t", header = T)
bl_b1d1=cbind(bl[,(1:9)], bl[,c(10,13)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b1d2=cbind(bl[,(1:9)], bl[,c(11,12)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b2d1=cbind(bl[,(1:9)], bl[,c(14,17)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b2d2=cbind(bl[,(1:9)], bl[,c(15,16)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b3d1=cbind(bl[,(1:9)], bl[,c(18,21)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b3d2=cbind(bl[,(1:9)], bl[,c(19,20)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b4d1=cbind(bl[,(1:9)], bl[,c(22,23)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b4d1=bl_b4d1[,-11]
bl_b4d2=cbind(bl[,(1:9)], bl[,c(22,23)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b4d2=bl_b4d2[,-12]
bl_b5d1=cbind(bl[,(1:9)], bl[,c(24,27)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b5d2=cbind(bl[,(1:9)], bl[,c(25,26)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b6d1=cbind(bl[,(1:9)], bl[,c(28,31)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b6d2=cbind(bl[,(1:9)], bl[,c(29,30)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])


write.table(EQP_b1d1, file='EQP_b1d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b1d2, file='EQP_b1d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b2d1, file='EQP_b2d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b2d2, file='EQP_b2d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b3d1, file='EQP_b3d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b3d2, file='EQP_b3d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b4d1, file='EQP_b4d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b4d2, file='EQP_b4d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b5d1, file='EQP_b5d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b5d2, file='EQP_b5d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b6d1, file='EQP_b6d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b6d2, file='EQP_b6d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b1d1, file='bl_b1d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b1d2, file='bl_b1d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b2d1, file='bl_b2d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b2d2, file='bl_b2d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b3d1, file='bl_b3d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b3d2, file='bl_b3d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b4d1, file='bl_b4d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b4d2, file='bl_b4d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b5d1, file='bl_b5d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b5d2, file='bl_b5d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b6d1, file='bl_b6d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b6d2, file='bl_b6d2.tsv', quote=FALSE, sep='\t')

add_blank_matches('EQP_b1d1.tsv', 'bl_b1d1.tsv')
add_blank_matches('EQP_b1d2.tsv', 'bl_b1d2.tsv')
add_blank_matches('EQP_b2d1.tsv', 'bl_b2d1.tsv')
add_blank_matches('EQP_b2d2.tsv', 'bl_b2d2.tsv')
add_blank_matches('EQP_b3d1.tsv', 'bl_b3d1.tsv')
add_blank_matches('EQP_b3d2.tsv', 'bl_b3d2.tsv')
add_blank_matches('EQP_b4d1.tsv', 'bl_b4d1.tsv')
add_blank_matches('EQP_b4d2.tsv', 'bl_b4d2.tsv')
add_blank_matches('EQP_b5d1.tsv', 'bl_b5d1.tsv')
add_blank_matches('EQP_b5d2.tsv', 'bl_b5d2.tsv')
add_blank_matches('EQP_b6d1.tsv', 'bl_b6d1.tsv')
add_blank_matches('EQP_b6d2.tsv', 'bl_b6d2.tsv')

EQP_b1d1_pbmatched=read.table("EQP_b1d1_pbmatched.tsv", sep="\t", header = T)
EQP_b1d2_pbmatched=read.table("EQP_b1d2_pbmatched.tsv", sep="\t", header = T)
EQP_b2d1_pbmatched=read.table("EQP_b2d1_pbmatched.tsv", sep="\t", header = T)
EQP_b2d2_pbmatched=read.table("EQP_b2d2_pbmatched.tsv", sep="\t", header = T)
EQP_b3d1_pbmatched=read.table("EQP_b3d1_pbmatched.tsv", sep="\t", header = T)
EQP_b3d2_pbmatched=read.table("EQP_b3d2_pbmatched.tsv", sep="\t", header = T)
EQP_b4d1_pbmatched=read.table("EQP_b4d1_pbmatched.tsv", sep="\t", header = T)
EQP_b4d2_pbmatched=read.table("EQP_b4d2_pbmatched.tsv", sep="\t", header = T)
EQP_b5d1_pbmatched=read.table("EQP_b5d1_pbmatched.tsv", sep="\t", header = T)
EQP_b5d2_pbmatched=read.table("EQP_b5d2_pbmatched.tsv", sep="\t", header = T)
EQP_b6d1_pbmatched=read.table("EQP_b6d1_pbmatched.tsv", sep="\t", header = T)
EQP_b6d2_pbmatched=read.table("EQP_b6d2_pbmatched.tsv", sep="\t", header = T)

EQP_matches=cbind(EQP_b1d1_pbmatched, EQP_b1d2_pbmatched, EQP_b2d1_pbmatched, EQP_b2d2_pbmatched, EQP_b3d1_pbmatched, EQP_b3d2_pbmatched, EQP_b4d1_pbmatched, EQP_b4d2_pbmatched, EQP_b5d1_pbmatched, EQP_b5d2_pbmatched, EQP_b6d1_pbmatched, EQP_b6d2_pbmatched)
colnames(EQP_matches)=make.unique(colnames(EQP_matches))
EQP_matches=cbind(EQP_matches$sample_min0, EQP_matches$blank_matches, EQP_matches$bf, EQP_matches$bf2, 
                  EQP_matches$sample_min0.1, EQP_matches$blank_matches.1, EQP_matches$bf.1, EQP_matches$bf2.1,
                  EQP_matches$sample_min0.2, EQP_matches$blank_matches.2, EQP_matches$bf.2, EQP_matches$bf2.2,
                  EQP_matches$sample_min0.3, EQP_matches$blank_matches.3, EQP_matches$bf.3, EQP_matches$bf2.3,
                  EQP_matches$sample_min0.3, EQP_matches$blank_matches.3, EQP_matches$bf.3, EQP_matches$bf2.3,
                  EQP_matches$sample_min0.4, EQP_matches$blank_matches.4, EQP_matches$bf.4, EQP_matches$bf2.4,
                  EQP_matches$sample_min0.5, EQP_matches$blank_matches.5, EQP_matches$bf.5, EQP_matches$bf2.5,
                  EQP_matches$sample_min0.6, EQP_matches$blank_matches.6, EQP_matches$bf.6, EQP_matches$bf2.6,
                  EQP_matches$sample_min0.7, EQP_matches$blank_matches.7, EQP_matches$bf.7, EQP_matches$bf2.7,
                  EQP_matches$sample_min0.8, EQP_matches$blank_matches.8, EQP_matches$bf.8, EQP_matches$bf2.8,
                  EQP_matches$sample_min0.9, EQP_matches$blank_matches.9, EQP_matches$bf.9, EQP_matches$bf2.9,
                  EQP_matches$sample_min0.10, EQP_matches$blank_matches.10, EQP_matches$bf.10, EQP_matches$bf2.10,
                  EQP_matches$sample_min0.11, EQP_matches$blank_matches.11, EQP_matches$bf.11, EQP_matches$bf2.11)
jolly=colnames(EQP_matches)=c("sample_min0", "blank_matches", "bf", "bf2", 
                              "sample_min0.1", "blank_matches.1", "bf.1", "bf2.1",
                              "sample_min0.2", "blank_matches.2", "bf.2", "bf2.2",
                              "sample_min0.3", "blank_matches.3", "bf.3", "bf2.3",
                              "sample_min0.3", "blank_matches.3", "bf.3", "bf2.3",
                              "sample_min0.4", "blank_matches.4", "bf.4", "bf2.4",
                              "sample_min0.5", "blank_matches.5", "bf.5", "bf2.5",
                              "sample_min0.6", "blank_matches.6", "bf.6", "bf2.6",
                              "sample_min0.7", "blank_matches.7", "bf.7", "bf2.7",
                              "sample_min0.8", "blank_matches.8", "bf.8", "bf2.8",
                              "sample_min0.9", "blank_matches.9", "bf.9", "bf2.9",
                              "sample_min0.10", "blank_matches.10", "bf.10", "bf2.10",
                              "sample_min0.11", "blank_matches.11", "bf.11", "bf2.11")
EQP_matches=data.frame(EQP_matches)
EQP_new <- filter_relevant(read.table('EQP.tsv',sep = "\t",header=TRUE))
jelly=EQP_new[-which((EQP_matches$sample_min0==0 & EQP_matches$blank_matches>=1) | (EQP_matches$blank_matches==1 & EQP_matches$bf<0) | (EQP_matches$bf<0 & EQP_matches$bf2<0) |
                       (EQP_matches$sample_min0.1==0 & EQP_matches$blank_matches.1>=1) | (EQP_matches$blank_matches.1==1 & EQP_matches$bf.1<0) | (EQP_matches$bf.1<0 & EQP_matches$bf2.1<0) |
                       (EQP_matches$sample_min0.2==0 & EQP_matches$blank_matches.2>=1) | (EQP_matches$blank_matches.2==1 & EQP_matches$bf.2<0) | (EQP_matches$bf.2<0 & EQP_matches$bf2.2<0) |
                       (EQP_matches$sample_min0.3==0 & EQP_matches$blank_matches.3>=1) | (EQP_matches$blank_matches.3==1 & EQP_matches$bf.3<0) | (EQP_matches$bf.3<0 & EQP_matches$bf2.3<0) |
                       (EQP_matches$sample_min0.4==0 & EQP_matches$blank_matches.4>=1) | (EQP_matches$blank_matches.4==1 & EQP_matches$bf.4<0) | (EQP_matches$bf.4<0 & EQP_matches$bf2.4<0) |
                       (EQP_matches$sample_min0.5==0 & EQP_matches$blank_matches.5>=1) | (EQP_matches$blank_matches.5==1 & EQP_matches$bf.5<0) | (EQP_matches$bf.5<0 & EQP_matches$bf2.5<0) |
                       (EQP_matches$sample_min0.6==0 & EQP_matches$blank_matches.6>=1) | (EQP_matches$blank_matches.6==1 & EQP_matches$bf.6<0) | (EQP_matches$bf.6<0 & EQP_matches$bf2.6<0) |
                       (EQP_matches$sample_min0.7==0 & EQP_matches$blank_matches.7>=1) | (EQP_matches$blank_matches.7==1 & EQP_matches$bf.7<0) | (EQP_matches$bf.7<0 & EQP_matches$bf2.7<0) |
                       (EQP_matches$sample_min0.8==0 & EQP_matches$blank_matches.8>=1) | (EQP_matches$blank_matches.8==1 & EQP_matches$bf.8<0) | (EQP_matches$bf.8<0 & EQP_matches$bf2.8<0) |
                       (EQP_matches$sample_min0.9==0 & EQP_matches$blank_matches.9>=1) | (EQP_matches$blank_matches.9==1 & EQP_matches$bf.9<0) | (EQP_matches$bf.9<0 & EQP_matches$bf2.9<0) |
                       (EQP_matches$sample_min0.10==0 & EQP_matches$blank_matches.10>=1) | (EQP_matches$blank_matches.10==1 & EQP_matches$bf.10<0) | (EQP_matches$bf.10<0 & EQP_matches$bf2.10<0) |
                       (EQP_matches$sample_min0.11==0 & EQP_matches$blank_matches.11>=1) | (EQP_matches$blank_matches.11==1 & EQP_matches$bf.11<0) | (EQP_matches$bf.11<0 & EQP_matches$bf2.11<0)),]

jelly=rbind(jelly,EQP_new[545,], EQP_new[1531,], EQP_new[1881,])
jel=t(jelly)
colnames(jel)=jel[1,]
tzel=jel[-c(1:9, 707:709),]

library(readr)
EQP_cl=read.delim(file="EQP_pb_cleaned_withMetaData.txt",sep = "\t",header=TRUE, dec = ".")
EQP_cle=t(EQP_cl)
colnames(EQP_cle)=EQP_cle[1,]
EQP_cle=EQP_cle[-1,]
EQP_cle=as.data.frame(EQP_cle)
EQP_cle$TrueNo=as.numeric(as.character(EQP_cle$TrueNo))
EQPc=EQP_cle[order(EQP_cle$TrueNo),]

tzelis=tzel[which(row.names(tzel) %in% row.names(EQPc)),]

num=function(x){
  which(row.names(tzelis)==row.names(EQPc)[x])
}
num_all=lapply(1:as.numeric(nrow(EQPc)), num)
numall=unlist(num_all)
tzelis=tzelis[numall,]

fab=function(x){
  batchX=EQPc[which(EQPc$BATCH==x),]
  #b1=batchX[,21:ncol(batchX)]
  b1=tzelis[row.names(tzelis) %in% row.names(batchX),]
  
  for (i in 1:ncol(b1)) {
    mi=b1[,i]
    mi=as.numeric(as.character(mi))
    for (j in 1:length(which(batchX$`F/B`=="forward"))) {
      if (mi[j]==0 | mi[length(mi)-j+1]==0){
        mi[j]=mi[length(mi)-j+1]=0
      } else {
        mi[j]=mean(c(mi[j], mi[length(mi)-j+1]))
      }}
    b1[,i]=mi
  }
  
  return(b1)
}

EQP_av=lapply(1:6, fab)
EQP_av_m=rbind(EQP_av[[1]],EQP_av[[2]],EQP_av[[3]],EQP_av[[4]],EQP_av[[5]],EQP_av[[6]])

make.num=function(x) {
  tzelis[,x]=as.numeric(as.character(tzelis[,x]))
}
tzelis=sapply(1:ncol(tzelis), make.num)
#replace wrongly corrected end of batch QCs with their true values
EQP_av_m[c(56,57,173,174,294,295,406,407,518,519,635,636),]=tzelis[c(56,57,173,174,294,295,406,407,518,519,635,636),1:ncol(tzelis)]
EQP_av_m=cbind(EQPc[,(1:20)], EQP_av_m)
EQP_av_m_m=EQP_av_m[-which(EQP_av_m$`F/B`=="backwards"),]

test=jel[which(row.names(EQP_av_m_m) %in% row.names(jel)),]
EQP_bp=t(EQP_av_m_m) #REMEMBER! <jel> has your variable names!

jack=t(EQP_bp)
jack=as.data.frame(jack)
pearl=jack[,(21:ncol(jack))]
indx <- sapply(pearl, is.factor)
pearl[indx] <- lapply(pearl[indx], function(x) as.numeric(as.character(x)))
pearl[pearl==0]=NA
#remove peaks not present in at least 80% of samples in each group
require("plyr")
count(jack, c("RUN"))
count(jack, c("group"))
#captain=pearl[,-which(colSums(is.na(pearl[which(jack$group=="QC plasmamonster"),]))>=3)]
captain2=pearl[,-which(colSums(is.na(pearl[which(jack$RUN=="run1"),]))>=11 | colSums(is.na(pearl[which(jack$RUN=="run2"),]))>=10 | colSums(is.na(pearl[which(jack$RUN=="run3"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run4"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run5"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run6"),]))>=11)]
captain=captain2[,-which(colSums(is.na(captain2[which(jack$group=="QC plasmamonster"),]))>=4)]
#---------START of adding IS------------
#captain=cbind(captain, pearl[,c(as.numeric(ncol(pearl)), as.numeric(ncol(pearl))-1, as.numeric(ncol(pearl))-2)])
#captain=captain[,-as.numeric(ncol(captain))]
#---------END of adding IS------------
#which(colSums(is.na(jack[which(jack$group=="CS"),]))>=5 & colSums(is.na(jack[which(jack$group=="PA"),]))>=13 & colSums(is.na(jack[which(jack$group=="PHT"),]))>=13 & colSums(is.na(jack[which(jack$group=="PPGL"),]))>=12)

test2=test[,-which(colSums(is.na(pearl[which(jack$RUN=="run1"),]))>=11 | colSums(is.na(pearl[which(jack$RUN=="run2"),]))>=10 | colSums(is.na(pearl[which(jack$RUN=="run3"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run4"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run5"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run6"),]))>=11)]
test2=test2[,-which(colSums(is.na(captain2[which(jack$group=="QC plasmamonster"),]))>=4)]
#test2=test2[,-as.numeric(ncol(test2))]
#test2=cbind(test2, test[,c(as.numeric(ncol(pearl)), as.numeric(ncol(pearl))-1, as.numeric(ncol(pearl))-2)])
#captain=captain[,-which(colSums(is.na(pearl[which(jack$group=="CS"),]))>=10 & colSums(is.na(pearl[which(jack$group=="PA"),]))>=26 & colSums(is.na(pearl[which(jack$group=="PPGL"),]))>=25 & colSums(is.na(pearl[which(jack$group=="PHT"),]))>=26)]
# calculate MV per sample, remove samples with >20% of total number of peaks
MV=rowSums(is.na(captain))
MV=MV/ncol(captain)
#captain=captain[-as.numeric(which(MV>0.3)),]
sparrow=jack#[-as.numeric(which(MV>0.3)),]
#TINT
captain[is.na(captain)]=0
sumcap=function(x){
  sum(captain[x,])
}
sums=sapply(1:as.numeric(nrow(captain)), sumcap)
captain[,as.numeric(ncol(captain))+1]=sums
tincap=captain/captain[,as.numeric(ncol(captain))]
tincap$V6343=NULL
captain=tincap
# PQN using HVs and using only peaks with no MV for calculating quotients
captain[captain==0]=NA
dauntless=captain
black=dauntless[which(sparrow$group=="HV" & sparrow$center=="Dresden (2)"),]
check1=apply(black, 2, median, na.rm=T)
check2=t(t(dauntless)/check1)
check3=apply(check2, 1, median, na.rm=T)
checkmate=captain/check3
checkmate[is.na(checkmate)]=0
#norrington=checkmate
#newcheck=checkmate*100000000
#checkmate=as.matrix(checkmate)
#test=sprintf("%.20f", checkmate)
#write.csv(test, "test.csv")
#write.csv(checkmate2, "test4.csv")
#test <- read_csv("Z:/NickB/all_neg2/BLANKS/results_XCMS_offline/test.csv")
#test4=as.matrix(test3)
#write.csv(newt, "test4.csv")
#test5 <- read_csv("Z:/NickB/all_neg2/BLANKS/results_XCMS_offline/test4.csv")
#Batch correction
checkmate.m=as.matrix(checkmate)
checkmate.m[which(checkmate.m==0)]=min(checkmate.m[which(checkmate.m>0)])
AV_i=checkmate.m
bn=c(rep(1,57), rep(2,62), rep(3,61), rep(4,53), rep(5,61), rep(6,58))
bn=as.factor(bn)
c_cor1=sweep(AV_i[which(bn==1),], 2, colMeans(AV_i[which(bn==1),]) ,'/')
c_cor2=sweep(AV_i[which(bn==2),], 2, colMeans(AV_i[which(bn==2),]) ,'/')
c_cor3=sweep(AV_i[which(bn==3),], 2, colMeans(AV_i[which(bn==3),]) ,'/')
c_cor4=sweep(AV_i[which(bn==4),], 2, colMeans(AV_i[which(bn==4),]) ,'/')
c_cor5=sweep(AV_i[which(bn==5),], 2, colMeans(AV_i[which(bn==5),]) ,'/')
c_cor6=sweep(AV_i[which(bn==6),], 2, colMeans(AV_i[which(bn==6),]) ,'/')

batcor=rbind(c_cor1, c_cor2, c_cor3, c_cor4, c_cor5, c_cor6)
batcor[checkmate==0]=0

checkmate=batcor

#RSD filter:30% threshold
hay=function(x) {
  sd(checkmate[which(sparrow$group=="QC plasmamonster"),x], na.rm = T)
}
sdQC_plasmamonster=sapply((1:(as.numeric(ncol(checkmate)))), hay)
bay=function(x) {
  mean(checkmate[which(sparrow$group=="QC plasmamonster"),x], na.rm = T)
}
meanQC_plasmamonster=sapply((1:(as.numeric(ncol(checkmate)))), bay)
cvQC_plasmamonster=sdQC_plasmamonster/meanQC_plasmamonster
norrington=checkmate[,which(cvQC_plasmamonster<0.3)]

#kNN missing value estimation
library(impute)
ncheck=norrington
ncheck[ncheck==0]=NA
ncheckO.imputed <- impute.knn(as.matrix(ncheck), k=10, rowmax = 1, colmax = 1, maxp=nrow(ncheck))
commodore=ncheckO.imputed$data
#glog
library("LMGene")
library(Biobase)
library(tools)
library(readr)

ENP_vanilla_forGLOG=commodore
PLASMA_QC_PQN=ENP_vanilla_forGLOG[which(sparrow$group=="QC plasmamonster"),]
QCled_PQN=as.matrix(t(PLASMA_QC_PQN))
QCr=apply(QCled_PQN, 1,as.numeric)
QCr=t(QCr)
QCmonster=as.factor(c(rep(1:2, each=(as.numeric(nrow(PLASMA_QC_PQN))/2)), 1))
#QCmonster=as.factor(c(rep(1:2, each=(as.numeric(nrow(PLASMA_QC_PQN))/2))))
QCdose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_QC_PQN)))))
QCled_list=list("monster"=QCmonster, "dose"=QCdose)
QCled.eS=neweS(QCr, QCled_list)
tranpar <- tranest(QCled.eS)
tranpar

PLASMA_PQN=ENP_vanilla_forGLOG
led_PQN=as.matrix(t(PLASMA_PQN))
r=apply(led_PQN, 1,as.numeric)
r=t(r)
monster=as.factor(c(1:as.numeric(nrow(PLASMA_PQN))))
dose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_PQN)))))
led_list=list("monster"=monster, "dose"=dose)
led.eS=neweS(r, led_list)
trled.eS <- transeS(led.eS, tranpar$lambda, tranpar$alpha)
kostakis=exprs(trled.eS)
colnames(kostakis)=as.factor(colnames(kostakis))
colnames(kostakis)=colnames(led_PQN)
final=t(kostakis)
FEATURE_NAMES=test2[1:9,which(cvQC_plasmamonster<0.3)]

write.csv(sparrow[,1:20], file="EQP_FINAL_Metadata_neg_from_R.csv")
EQP_FINAL_Metadata_neg_from_R <- read_csv("D:/Office 2010 Toolkit and EZ-Activator v2.2.3/metabolomics hypertension/QTOF_EQP_ENSAT-HT/EQP_FINAL_Metadata_neg_from_R.csv")
EQP_METADATA_NEG=EQP_FINAL_Metadata_neg_from_R[order(EQP_FINAL_Metadata_neg_from_R$X1),]
EQP_NEG=final[order(row.names(final)),]
row.names(EQP_NEG)==EQP_METADATA_NEG$X1

library(mixOmics)
X=EQP_NEG
row.names(X)=EQP_METADATA_NEG$`sample name 2`
pca.EQP = pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = EQP_METADATA_NEG$GROUP, title = "NEG QTOF", ind.names = F,legend = T,cex=4, pch = "sphere", style = "3d")

X=EQP_NEG[-which(EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV"),]
pca.EQP = pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = EQP_METADATA_NEG$GROUP[-which(EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV")], title = "NEG QTOF", ind.names = F,legend = T,cex=4, pch = "sphere", style = "3d")
plotIndiv(pca.EQP, group = EQP_METADATA_NEG$`CENTER ID`[-which(EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV")], title = "NEG QTOF", ind.names = F,legend = T,cex=4, pch = "sphere", style = "3d")

WHYF=data.frame(EQP_METADATA_NEG[-which(EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV"),])
sa=as.numeric(WHYF$`SAMPLE.AGE`)
kruskal.test(sa~WHYF$CENTER.ID)
#p-value < 2.2e-16
sa[which(sa<median(sa))]="<median"
sa[-which(sa=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = sa, legend = T, pch=c(19), cex = 5, title="(c): PCA SAMPLE AGE")
cor.test(pI$df$x, sa, method = "kendall")
#p-value = 1.329e-11
cor.test(pI$df$x, sa, method = "spearman")
#p-value = 1.997e-12
cor.test(pI$df$y, sa, method = "kendall")
#p-value = 0.5918
cor.test(pI$df$y, sa, method = "spearman")
#p-value = 0.5878
SEQ=as.numeric(WHYF$QTOF.No)
kruskal.test(SEQ~WHYF$CENTER.ID)
#p-value = 0.03315
SEQ[which(SEQ<median(SEQ))]="<median"
SEQ[-which(SEQ=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = SEQ, legend = T, pch=c(19), cex = 5, title="(c): PCA RUN ORDER")
cor.test(pI$df$x, SEQ, method = "kendall")
#p-value = 0.9149
cor.test(pI$df$x, SEQ, method = "spearman")
#p-value = 0.8621
cor.test(pI$df$y, SEQ, method = "kendall")
#p-value = 0.4544
cor.test(pI$df$y, SEQ, method = "spearman")
#p-value = 0.445
RUN=as.numeric(WHYF$BATCH)
kruskal.test(RUN~WHYF$CENTER.ID)
#p-value = 0.6048
RUN[which(RUN<median(RUN))]="<median"
RUN[-which(RUN=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = RUN, legend = T, pch=c(19), cex = 5, title="(c): PCA BATCH")
cor.test(pI$df$x, RUN, method = "kendall")
#p-value =  0.7987
cor.test(pI$df$x, RUN, method = "spearman")
#p-value = 0.8128
cor.test(pI$df$y, RUN, method = "kendall")
#p-value = 0.9312
cor.test(pI$df$y, RUN, method = "spearman")
#p-value = 0.8796


X1=X
colnames(X1)=FEATURE_NAMES[2,]
pca.EQP = pca(X1, ncomp = 10, center = TRUE, scale = FALSE)
plotVar(pca.EQP, cex=3)

#final=cbind(sparrow[,1:20], final)

WHYF=data.frame(EQP_METADATA_NEG)
WHYF=WHYF[which(WHYF$GROUP=="PPGL" | WHYF$GROUP=="CS" | WHYF$GROUP=="PA" | WHYF$GROUP=="PHT"),]
ex=X
Y=WHYF$GROUP
length(which(Y=="CS"))
length(which(Y=="PA"))
length(which(Y=="PHT"))
length(which(Y=="PPGL"))

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="CS")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="CS"),])))
    folds2=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds3=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds4=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 7, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 7, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = F)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.train$choice.keepX)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    MCs=jack[-which(Prediction==Y),(1:14)]
    CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.CS=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.CS=lapply((1:n), each.CS)
all.CS=unlist(all.CS)

each.PA=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

each.PPGL=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)






plsda.select <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
system.time(perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = f, auc = F, progressBar = T, nrepeat = n, dist="mahalanobis.dist"))
#600sec
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
system.time(tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = f, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = n, scale = F, cpus = 6))
#80sec
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]
test.predict <- predict(splsda.select, X, dist = "mahalanobis.dist")
B.hat_allvall1=test.predict$B.hat[,,ncomp]
rich=B.hat_allvall1#[,1]
rich=rich[which(rowSums(rich)>0 | rowSums(rich)<0),]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)



allvall1=list(plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)
allvall2=list(allconf, all.BER, all.CS, all.PA , all.PPGL , all.PHT)


hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    fit = glmnet(X, Y, family = "multinomial", alpha=1, standardize = F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 7, type.measure = "class", standardize = F)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(cvfit1, X, type = "class", s = "lambda.min")     #change distance
    confusion.mat = get.confusion_matrix(truth = Y, predicted = test.predict)
    BER.res = get.BER(confusion.mat)
    sensitivity=confusion.mat[1,1]/sum(confusion.mat[1,])
    specificity=confusion.mat[2,2]/sum(confusion.mat[2,])
    rich=coef(cvfit1, s = "lambda.min")
    MCs=jack[-which(test.predict==Y),(1:14)]
    CCs=jack[which(test.predict==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, test.predict, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.CS=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.CS=lapply((1:n), each.CS)
all.CS=unlist(all.CS)

each.PA=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

each.PPGL=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)




#X=ex
#Y[which(Y=="AHT")]="1"
#Y[which(Y=="PHT")]="0"
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 8, type.measure = "class", standardize=FALSE)
richCS=coef(cvfit1, s = "lambda.min")$CS
richPA=coef(cvfit1, s = "lambda.min")$PA
richPHT=coef(cvfit1, s = "lambda.min")$PHT
richPPGL=coef(cvfit1, s = "lambda.min")$PPGL
rich=cbind(richCS, richPA, richPHT, richPPGL)
colnames(rich)=c("CS", "PA", "PHT", "PPGL")
rich=rich[-1,]
rich=rich[-which(rowSums(rich)==0),]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)

allvall3=list(cvfit1, rich)
allvall4=list(allconf, all.BER, all.CS, all.PA , all.PPGL , all.PHT)

length(which(rich[,1]==0 & rich[,2]==0 & rich[,3]==0 & rich[,4]==0))
length(which(rowSums(rich)==0))

sex=as.factor(WHYF$GENDER)
levels(sex)=c("0", "1")
levels(sex)=as.numeric(levels(sex))
sex=as.numeric(sex)
sex[sex==1]=0
sex[sex==2]=1
age=as.numeric(WHYF$PATIENT.AGE)
ex=cbind(ex, sex, age)

X=ex
Y=WHYF$GROUP
#colnames(X)=round(as.numeric(colnames(X)),3)
pca.EQP=mixOmics::pca(X, scale = T, center = T, ncomp = 10)
plotIndiv(pca.EQP, group = Y, legend = T, style = "3d", pch="sphere")
plotIndiv(pca.EQP, group = WHYF$CENTER.ID, legend = T, style = "3d", pch="sphere")
length(which(Y=="CS"))
length(which(Y=="PA"))
length(which(Y=="PHT"))
length(which(Y=="PPGL"))

sa=as.numeric(WHYF$`SAMPLE.AGE`)
kruskal.test(sa~WHYF$CENTER.ID)
#p-value < 2.2e-16
sa[which(sa<median(sa))]="<median"
sa[-which(sa=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = sa, legend = T, pch=c(19), cex = 5, title="(c): PCA SAMPLE AGE")
cor.test(pI$df$x, sa, method = "kendall")
#p-value = 8.157e-05
cor.test(pI$df$x, sa, method = "spearman")
#p-value = 7.063e-05
cor.test(pI$df$y, sa, method = "kendall")
#p-value = 4.133e-06
cor.test(pI$df$y, sa, method = "spearman")
#p-value = 7.584e-07
SEQ=as.numeric(WHYF$QTOF.No)
kruskal.test(SEQ~WHYF$CENTER.ID)
#p-value = 0.03315
SEQ[which(SEQ<median(SEQ))]="<median"
SEQ[-which(SEQ=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = SEQ, legend = T, pch=c(19), cex = 5, title="(c): PCA RUN ORDER")
cor.test(pI$df$x, SEQ, method = "kendall")
#p-value = 0.1179
cor.test(pI$df$x, SEQ, method = "spearman")
#p-value = 0.1197
cor.test(pI$df$y, SEQ, method = "kendall")
#p-value = 0.2799
cor.test(pI$df$y, SEQ, method = "spearman")
#p-value = 0.2713
RUN=as.numeric(WHYF$BATCH)
kruskal.test(RUN~WHYF$CENTER.ID)
#p-value = 0.6048
RUN[which(RUN<median(RUN))]="<median"
RUN[-which(RUN=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = RUN, legend = T, pch=c(19), cex = 5, title="(c): PCA BATCH")
cor.test(pI$df$x, RUN, method = "kendall")
#p-value = 0.6953
cor.test(pI$df$x, RUN, method = "spearman")
#p-value = 0.6722
cor.test(pI$df$y, RUN, method = "kendall")
#p-value = 0.7117
cor.test(pI$df$y, RUN, method = "spearman")
#p-value = 0.7038

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="CS")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="CS"),])))
    folds2=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds3=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds4=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    plsda.train <- plsda(X, Y, ncomp = 10, scale=T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 7, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 7, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = T)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = T, keepX = tune.splsda.train$choice.keepX)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    MCs=jack[-which(Prediction==Y),(1:14)]
    CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.CS=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.CS=lapply((1:n), each.CS)
all.CS=unlist(all.CS)

each.PA=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

each.PPGL=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)







X=ex
plsda.select <- plsda(X, Y, ncomp = 10, scale=T)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = 8, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = 8, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 50, scale = T)
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = T, keepX = tune.splsda.select$choice.keepX)
test.predict <- predict(splsda.select, X, dist = "mahalanobis.dist")
B.hat_allvall5=test.predict$B.hat[,,ncomp]
rich=B.hat_allvall5#[,1]
rich=rich[which(rowSums(rich)>0 | rowSums(rich)<0),]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)



allvall5=list(plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)
allvall6=list(allconf, all.BER, all.CS, all.PA , all.PPGL , all.PHT)


hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    fit = glmnet(X, Y, family = "multinomial", alpha=1, standardize = T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 7, type.measure = "class", standardize = T)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(cvfit1, X, type = "class", s = "lambda.min")     #change distance
    confusion.mat = get.confusion_matrix(truth = Y, predicted = test.predict)
    BER.res = get.BER(confusion.mat)
    sensitivity=confusion.mat[1,1]/sum(confusion.mat[1,])
    specificity=confusion.mat[2,2]/sum(confusion.mat[2,])
    rich=coef(cvfit1, s = "lambda.min")
    MCs=jack[-which(test.predict==Y),(1:14)]
    CCs=jack[which(test.predict==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, test.predict, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.CS=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.CS=lapply((1:n), each.CS)
all.CS=unlist(all.CS)

each.PA=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

each.PPGL=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)




X=as.matrix(ex)
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 8, type.measure = "class", standardize=T)
richCS=coef(cvfit1, s = "lambda.min")$CS
richPA=coef(cvfit1, s = "lambda.min")$PA
richPHT=coef(cvfit1, s = "lambda.min")$PHT
richPPGL=coef(cvfit1, s = "lambda.min")$PPGL
rich=cbind(richCS, richPA, richPHT, richPPGL)
colnames(rich)=c("CS", "PA", "PHT", "PPGL")
rich=rich[-1,]
rich=rich[-which(rowSums(rich)==0),]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)

allvall7=list(cvfit1, rich)
allvall8=list(allconf, all.BER, all.CS, all.PA , all.PPGL , all.PHT)
